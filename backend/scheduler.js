const cron = require("node-cron");
const { PythonShell } = require('python-shell');
const path = require('path');

// May need specific command to set it up?

function runAPIJob() {
    return new Promise((resolve, reject) => {
        const scriptName = 'API_JOB.py';
        // API_JOB.py is in the root of senior-proj, scheduler.js is in backend.
        const scriptDir = path.resolve(__dirname, '..'); // Resolves to /home/coleoliva/senior-proj
        const fullScriptPath = path.join(scriptDir, scriptName);

        let options = {
            mode: 'text',
            pythonPath: 'python3', // Ensure this is the correct path to your python3 executable
            pythonOptions: ['-u'], // Unbuffered stdout/stderr
            scriptPath: scriptDir, // Directory containing the Python script
            // args: ['value1', 'value2'] // If your script needs arguments
        };

        console.log(`Attempting to run Python script: ${scriptName} from directory ${scriptDir}`);
        let pyshell = new PythonShell(scriptName, options);

        let outputMessages = [];
        let errorMessages = [];

        pyshell.on('message', (message) => {
            console.log(`Python script output: ${message}\n`);
            outputMessages.push(message);
        });

        pyshell.on('stderr', (stderr) => {
            console.error(`Python script stderr: ${stderr}`);
            errorMessages.push(stderr);
        });

        pyshell.once('error', (err) => {
            console.error(`PythonShell failed to start or crashed: ${err}`);
            // Ensure reject is called only once if not already terminated
            if (!pyshell.terminated) {
                reject(new Error(`PythonShell execution error: ${err.message || err}`));
            }
        });

        pyshell.end((err, code, signal) => {
            if (err) {
                console.error(`PythonShell process error: ${err}`);
                reject(new Error(`Python script execution failed: ${err.message || err}. Output: ${outputMessages.join('\n')}. Stderr: ${errorMessages.join('\n')}`));
                return;
            }
            // A non-zero exit code from the Python script often indicates an error within the script
            if (code !== 0) {
                console.warn(`Python script finished with non-zero exit code: ${code}, signal: ${signal}`);
                reject(new Error(`Python script exited with code ${code}. Output: ${outputMessages.join('\n')}. Stderr: ${errorMessages.join('\n')}`));
                return;
            }
            console.log(`Python script finished successfully. Code: ${code}, Signal: ${signal}`);
            resolve({ success: true, message: "Python script executed successfully.", output: outputMessages, stderr: errorMessages });
        });
    });
}

// Schedule the task to run every Monday at 2 AM
cron.schedule('0 2 * * 1', () => {
    console.log('Running scheduled API job...');
    runAPIJob()
        .then(result => console.log('Scheduled job completed successfully:', result.message))
        .catch(error => console.error('Scheduled job failed:', error.message));
});

console.log('Scheduler setup. API Job will run on schedule.');

module.exports = { runAPIJob }; // Export the function