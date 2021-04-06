import { app, BrowserWindow, Menu, shell } from 'electron';

// Handle creating/removing shortcuts on Windows when installing/uninstalling.
if (require('electron-squirrel-startup')) { // eslint-disable-line global-require
  app.quit();
}
const path = require('path')

// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
let mainWindow;
let treeWindow;

let params; // User parameters set in main window

let template = [{
  label: 'File',
  submenu: [{
    label: 'Load options',
    accelerator: 'CmdOrCtrl+L',
    click: (item, focusedWindow) => {
      focusedWindow.send('load-options');
    }
  }, {
    label: 'Save options',
    accelerator: 'CmdOrCtrl+S',
    click: (item, focusedWindow) => {
      focusedWindow.send('save-options');
    }
  }, {
    label: 'Reset default option',
    accelerator: 'CmdOrCtrl+R',
    click: (item, focusedWindow) => {
      focusedWindow.send('reset-options');
    }
  }, 
  {
    type: 'separator'
  }, {
    label: 'Exit',
    accelerator: 'CmdOrCtrl+Q',
    role: 'quit'
  }]
}, {
  label: 'View',
  submenu: [{
    label: 'Toggle Full Screen',
    accelerator: (() => {
      if (process.platform === 'darwin') {
        return 'Ctrl+Command+F'
      } else {
        return 'F11'
      }
    })(),
    click: (item, focusedWindow) => {
      if (focusedWindow) {
        focusedWindow.setFullScreen(!focusedWindow.isFullScreen())
      }
    }
  }]
}, {
  label: 'Help',
  role: 'help',
  submenu: [{
    label: 'Getting started',
    click: () => {
      shell.openExternal('https://github.com/524D/compareMS2')
    }
  }]
}]

const menu = Menu.buildFromTemplate(template)
Menu.setApplicationMenu(menu)


const createWindow = () => {
  // Create the browser window.
  mainWindow = new BrowserWindow({
    width: 700,
    height: 650,
    webPreferences: {
      nodeIntegration: true,
      enableRemoteModule: true,
      contextIsolation: false,  // without this, we can't open new windows
      preload: path.join(__dirname, 'preload.js')
    }

  });

  // and load the index.html of the app.
  mainWindow.loadURL(`file://${__dirname}/index.html`);

  if (typeof process.env.CMPMS2_DEBUG !== 'undefined') {
      // Open the DevTools.
      mainWindow.webContents.openDevTools();
  } 

  // Emitted when the window is closed.
  mainWindow.on('closed', () => {
    // Dereference the window object, usually you would store windows
    // in an array if your app supports multi windows, this is the time
    // when you should delete the corresponding element.
    mainWindow = null;
  });
};

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', createWindow);

// Quit when all windows are closed.
app.on('window-all-closed', () => {
  // On OS X it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  // On OS X it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (mainWindow === null) {
    createWindow();
  }
});

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and import them here.
const {ipcMain, dialog} = require('electron')

ipcMain.on('open-dir-dialog', (event) => {
  const files = dialog.showOpenDialogSync(mainWindow, {
    properties: ['openDirectory']
  });
  if (files) {
      mainWindow.send('selected-directory', files)
  }
})

ipcMain.on('open-speciesfile-dialog', (event) => {
  const files = dialog.showOpenDialogSync(mainWindow, {
    filters: {
      filters: [
        { name: 'Text file', extensions: ['txt'] },
        { name: 'All Files', extensions: ['*'] }
      ]
    },
    properties: ['openFile']
  });
  if (files) {
      mainWindow.send('selected-speciesfile', files)
  }
})

// Display tree windows and send params
ipcMain.on('maketree', (event, args) => {
  const modalPath = path.join('file://', __dirname, '/tree.html')
  params = args;
  treeWindow = new BrowserWindow({
    width: 800,
    height: 700,
    parent: mainWindow,
    modal: true,
    webPreferences: {
        nodeIntegration: true,
        enableRemoteModule: true,
        contextIsolation: false,  // without this, we can't open new windows
        preload: path.join(__dirname, 'preload.js')
    }
  })
  treeWindow.on('close', () => { treeWindow = null })
  treeWindow.removeMenu();
  treeWindow.loadURL(modalPath);
  if (typeof process.env.CMPMS2_DEBUG !== 'undefined') {
    // Open the DevTools.
    treeWindow.webContents.openDevTools();
  } 

  treeWindow.show();
})

// Send parameters to tree window
ipcMain.on('get-userparms', () => {
   treeWindow.send('userparams', params);
})
