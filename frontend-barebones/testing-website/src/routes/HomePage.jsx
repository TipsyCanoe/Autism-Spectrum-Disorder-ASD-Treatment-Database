// general code for home page

import React from "react";
import "../index.css";

const HomePage = () => {
  return (
    <div>
      {/* Header Section */}
      <header className="bg-blue-600 text-white py-4">
        <div className="container flex justify-between items-center">
          <h1 className="text-2xl font-bold">My Website</h1>
          <nav>
            <ul className="flex space-x-6">
              <li>
                <a href="#" className="hover:underline">
                  Home
                </a>
              </li>
              <li>
                <a href="#" className="hover:underline">
                  About
                </a>
              </li>
              <li>
                <a href="#" className="hover:underline">
                  Contact
                </a>
              </li>
            </ul>
          </nav>
        </div>
      </header>

      {/* Hero Section */}
      <section className="bg-gray-100 py-20">
        <div className="container text-center">
          <h2 className="text-4xl font-bold text-gray-800">
            Welcome to My Website
          </h2>
          <p className="text-gray-600 mt-4">
            A simple homepage layout to showcase the website's look and feel.
          </p>
          <button className="mt-6 px-6 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700">
            Get Started
          </button>
        </div>
      </section>

      {/* Features Section */}
      <section className="bg-white py-16">
        <div className="container">
          <h3 className="text-2xl font-bold text-gray-800 text-center mb-8">
            Features
          </h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h4 className="text-xl font-semibold">Feature 1</h4>
              <p className="text-gray-600 mt-2">
                Description of feature 1 goes here. Brief and to the point.
              </p>
            </div>
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h4 className="text-xl font-semibold">Feature 2</h4>
              <p className="text-gray-600 mt-2">
                Description of feature 2 goes here. Explain its benefit.
              </p>
            </div>
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h4 className="text-xl font-semibold">Feature 3</h4>
              <p className="text-gray-600 mt-2">
                Description of feature 3 goes here. Short and informative.
              </p>
            </div>
          </div>
        </div>
      </section>

      {/* Footer Section */}
      <footer className="bg-gray-800 text-white text-center py-4">
        <p>&copy; 2025 My Website. All rights reserved.</p>
      </footer>
    </div>
  );
};

export default HomePage;
