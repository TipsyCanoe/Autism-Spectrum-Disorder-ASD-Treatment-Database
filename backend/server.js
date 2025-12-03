// Server entry point for the job scheduler
console.log('Starting Job Scheduler Server...');

// Import the scheduler
require('./scheduler');

// Keep the process alive
setInterval(() => {
    // Heartbeat every hour
}, 3600000);
