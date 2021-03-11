const {ipcRenderer} = nodeRequire('electron');
const app = nodeRequire('electron').remote.app;
const path = nodeRequire('path');
const { spawn } = nodeRequire('child_process');
const lineReader = nodeRequire('line-reader');

let tree = d3.layout.phylotree()
  // create a tree layout object
  .svg(d3.select("#tree_display"));
  // render to this SVG element
let svgTree = document.getElementById('tree_display');
//tree.size([svgTree.clientWidth,svgTree.clientHeight]);
tree.size([800,600]);
tree.font_size(15);
tree.options({'left-right-spacing' : "fit-to-size",
              'top-bottom-spacing': "fit-to-size",
              'zoom': true,}, true);

$("#layout").on("click", function(e) {
    tree.radial($(this).prop("checked")).placenodes().update();
  });

function escapeHtml(unsafe) {
    return unsafe
         .replace(/&/g, "&amp;")
         .replace(/</g, "&lt;")
         .replace(/>/g, "&gt;")
         .replace(/"/g, "&quot;")
         .replace(/'/g, "&#039;");
}

let file1_idx;
let file2_idx;
let paramsGlobal;  // To save memory in recursive call, we store these in global variables
let mgfFilesGlobal;
let compareMS2exe;
let compToDistExe;
let compResultListFile;
const myPath = app.getAppPath();

if (navigator.platform=='Linux x86_64') {
    compareMS2exe = path.join(myPath, 'external_binaries', 'compareMS2');
    compToDistExe = path.join(myPath, 'external_binaries', 'compareMS2_to_distance_matrices');
} else if ( (navigator.platform=='Win64') || (navigator.platform=='Win32')) {
    compareMS2exe = path.join(myPath, 'external_binaries', 'compareMS2.exe');
    compToDistExe = path.join(myPath, 'external_binaries', 'compareMS2_to_distance_matrices.exe');
}
else {
    document.body.innerHTML = "<H1>This app runs only on 64 bit Windows or 64 bit Linux Intel/AMD</H1>";
}

function compareNext() {
    let act=document.getElementById('activity');

    if (file1_idx >= mgfFilesGlobal.length) {
        act.innerHTML = 'Finished';
    }
    else
    {
        act.innerHTML = 'Comparing<br/>' + escapeHtml(mgfFilesGlobal[file1_idx]) + '<br/>' + mgfFilesGlobal[file2_idx];
        let cmpFile = path.join(paramsGlobal.mgfDir, "cmp_"+file1_idx+"_"+file2_idx+".txt");
        const cmp_ms2 = spawn(compareMS2exe,
        ['-1', mgfFilesGlobal[file1_idx],
        '-2', mgfFilesGlobal[file2_idx],
        '-c', paramsGlobal.cutoff,
        '-p', paramsGlobal.precMassDiff,
        '-w', paramsGlobal.chromPeakW,
        '-o', cmpFile,
        ]);
    
        cmp_ms2.stdout.on('data', (data) => {
            data = escapeHtml(data.toString());
            data = data.replace(/(?:\r\n|\r|\n)/g, '<br>');
            data = data.replace(/(?: )/g, '&nbsp;');
            document.getElementById('stdout').innerHTML += data;
            });
            
        cmp_ms2.stderr.on('data', (data) => {
            console.error(`stderr: ${data}`);
        });
            
        cmp_ms2.on('error', (data) => {
            console.error('Error running compareMS2');
            act.innerHTML = 'Error running compareMS2';
        });
            
        cmp_ms2.stderr.on('exit', (code, signal) => {
            console.error('Error running compareMS2');
            act.innerHTML = 'Error running compareMS2';
        });
            
        cmp_ms2.on('close', (code) => {
            fs.appendFileSync(compResultListFile, cmpFile + "\n");
            file2_idx++;
            if (file2_idx<file1_idx) {
                // If row is not finished, schedule next comparison
                setTimeout(function() {compareNext();}, 0);
            }
            else {
                // Finished new row, create tree
                makeTree();
            }
        });
    }
}

function makeTree() {
    let act=document.getElementById('activity');
    act.innerHTML = 'Creating tree';

    let cmdArgs = ['-i', compResultListFile,
    '-o', path.join(paramsGlobal.mgfDir, paramsGlobal.outBasename) ,
    '-c', paramsGlobal.cutoff,
    '-m'  // Generate MEGA format
    ]
    let s2s = paramsGlobal.s2sFile;
    // If the file to species mapping file exists, use it
    if (fs.existsSync(s2s) && fs.lstatSync(s2s).isFile()) {
        cmdArgs.push('-x', s2s)
    }
    else {
    // FIXME: compareMS2_to_distance_matrices doesn't work without sample2species file,
    // so assume it is in the data dir if not specified
        cmdArgs.push('-x', path.join(paramsGlobal.mgfDir, 'sample_to_species.txt'));
    }
    //                const c2d = spawn('echo', cmdArgs);
    const c2d = spawn(compToDistExe, cmdArgs);
    c2d.stdout.on('data', (data) => {
        data = escapeHtml(data.toString());
        data = data.replace(/(?:\r\n|\r|\n)/g, '<br>');
        data = data.replace(/(?: )/g, '&nbsp;');
        document.getElementById('stdout').innerHTML += data;
        });
        
    c2d.stderr.on('data', (data) => {
        console.error(`stderr: ${data}`);
    });
        
    c2d.on('error', (data) => {
        console.error('Error running compareMS2_to_distance_matrices');
        act.innerHTML = 'Error running compareMS2_to_distance_matrices';
    });
        
    c2d.stderr.on('exit', (code, signal) => {
        console.error('Error running compareMS2_to_distance_matrices');
        act.innerHTML = 'Error running compareMS2_to_distance_matrices';
    });
        
    c2d.on('close', (code) => {
        act.innerHTML = 'Computing tree';
        // Extract matrix and names from compareMS2_to_distance_matrices output
        let parseState = 'init';
        const reSpecies = /^QC\s+(.+)\s+([0-9\.]+)$/;
        const reMatrix = /^[0-9. \t]+$/;
        const reMatrixCapt = /([0-9\.]+)/g;
        let labels = [];
        let matrix = []; // Will be filled with rows -> 2D matrix
        matrix[0] = []; // First element must be empty
        const df = path.join(paramsGlobal.mgfDir, paramsGlobal.outBasename+'_distance_matrix.meg');
        lineReader.eachLine(df, (line, last) => {
            // Reading line by line
            if ( (parseState == 'init') ||
                    (parseState == 'labels') ) {
                let s = line.match(reSpecies);
                if ( (s) && (s.length != 0) ){
                    parseState = 'labels';
                    labels.push(s[1]);
                    // TODO: use CQ value
                } else if (parseState == 'labels') {
                    parseState = 'matrix';
                }
            }
            if (parseState == 'matrix') {
                if (reMatrix.test(line)) {
                    let row = line.match(reMatrixCapt);
                    // First one contains whole string, remove
                    // row.shift();
                    // Convert strings to numbers
                    row = row.map(x=>+x)
                    matrix.push(row);
                }
            }
            // Create new tree when files has finished loading
            if (last) {
                // Convert matrix and names into Newick format
                act.innerHTML = 'Showing tree';
                let newick = UPGMA(matrix, labels);
                console.log('newick', newick);
                d3.select("#tree_display").selectAll("*").remove();
                tree(newick)
                    .layout();
                file2_idx=0;

                file1_idx++;
                document.getElementById('stdout').innerHTML = '';
                // Start next comparison (if any)
                setTimeout(function() {compareNext();}, 0);
            }
        });
    });
}

function runCompare(params) {
    // TODO: sanitize params
    mgfFilesGlobal = getMgfFiles(params.mgfDir);
    // compareMS2 executables need local filenames, so change default dir
    process.chdir(params.mgfDir);
    // TODO: Sort files according to setting
    paramsGlobal = params;
    file1_idx = 1;
    file2_idx = 0;

    // Create empty comparison list file
    compResultListFile = path.join(paramsGlobal.mgfDir,'cmp_list.txt');
    fs.closeSync(fs.openSync(compResultListFile, 'w'))
    compareNext();
}

// Receive parameters set in the main window
ipcRenderer.on('userparams', (event, params) => {
    runCompare(params);
})

// Notify main process that we are ready to receive parameters
ipcRenderer.send('get-userparms');

