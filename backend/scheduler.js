const cron = require('node-cron');
const { exec } = require('child_process');
const path = require('path');

console.log('Scheduler service started');

// Schedule task to run every Sunday at 7:00 AM
// Cron format: Minute Hour DayOfMonth Month DayOfWeek
// 0 7 * * 0 = At 07:00 on Sunday
cron.schedule('0 7 * * 0', () => {
    console.log('Running weekly data pull...');
    
    // Path to API_JOB.py in the root directory
    const scriptPath = path.join(__dirname, '../API_JOB.py');
    
    // Execute the python script
    exec(`python3 ${scriptPath}`, {
        cwd: path.join(__dirname, '..') // Run from root directory
    }, (error, stdout, stderr) => {
        if (error) {
            console.error(`Error executing API_JOB.py: ${error}`);
            return;
        }
        if (stderr) {
            console.error(`API_JOB.py stderr: ${stderr}`);
        }
        console.log(`API_JOB.py output: ${stdout}`);
    });
});

console.log('Job scheduled for Sunday at 7:00 AM.');
