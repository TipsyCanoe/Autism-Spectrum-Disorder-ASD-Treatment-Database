const cron = require("node-cron");
const {PythonShell} = require('python-shell') 
const path = require('path');

// May need specific command to set it up?

function runAPIJob() {
    const scriptPath = path.join(__dirname, 'API_JOB.py');

    let options = {
        pythonPath: 'python3',
        scriptPath: __dirname,  // Automatically uses the current script directory
    };

    let pyshell = new PythonShell(scriptPath, options);

    pyshell.on('message', (message) => {
        console.log(`Output: ${message}`);
    });

    pyshell.on('error', (error) => {
        console.error(`Error: ${error}`);
    });

    pyshell.end((err, code, signal) => {
        if (err) console.error(`Error: ${err}`);
        console.log(`Python script finished with code: ${code}, signal: ${signal}`);
    });
}

// Schedule the task to run every Monday at 2 AM
cron.schedule('0 2 * * 1', runAPIJob);

console.log('API Job running...');