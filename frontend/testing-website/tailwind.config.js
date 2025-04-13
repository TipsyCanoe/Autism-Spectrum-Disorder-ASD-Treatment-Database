/** @type {import('tailwindcss').Config} */
module.exports = {
  content: [
    "./src/**/*.{js,jsx,ts,tsx}", // Ensure this covers your project structure
  ],
  theme: {
    extend: {
      colors: {
        "navbar-blue": "#2c3e50", // Custom name for the background
        "link-hover-blue": "#3498db", // Custom name for the hover background
        custom: "#3498db", // Assuming 'bg-custom' from HomePage uses this
      },
      spacing: {
        3.75: "0.9375rem", // Custom spacing for 15px (15 / 16 = 0.9375)
        1.25: "0.3125rem", // Custom spacing for 5px (5 / 16 = 0.3125) - if needed for radius
      },
      borderRadius: {
        "md-5": "5px", // Custom border radius if needed instead of md
      },
    },
  },
  plugins: [],
};
