import React from "react";
import "../index.css";

const HomePage = () => {
  return (
    <div>
      {/* Hero Section */}
      <section className="bg-gray-100 py-20">
        <div className="container mx-auto text-center">
          <h2 className="text-4xl font-bold text-gray-800">
            Welcome to the Autism Resources Database
          </h2>
          <p className="text-gray-600 mt-4">
            A comprehensive database of resources for professionals and
            families.
          </p>
          <button className="mt-6 px-6 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700">
            Explore Resources
          </button>
        </div>
      </section>

      {/* Features Section */}
      <section className="bg-white py-16">
        <div className="container mx-auto">
          <h3 className="text-2xl font-bold text-gray-800 text-center mb-8">
            Features
          </h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h4 className="text-xl font-semibold">Extensive Database</h4>
              <p className="text-gray-600 mt-2">
                Access a wide range of resources tailored for autism support.
              </p>
            </div>
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h4 className="text-xl font-semibold">Professional Guidance</h4>
              <p className="text-gray-600 mt-2">
                Find expert advice and best practices for autism care.
              </p>
            </div>
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h4 className="text-xl font-semibold">Family Support</h4>
              <p className="text-gray-600 mt-2">
                Discover resources and support networks for families.
              </p>
            </div>
          </div>
        </div>
      </section>

      {/* Footer Section */}
      <footer className="bg-gray-800 text-white text-center py-4">
        <p>&copy; 2025 Autism Resources Database. All rights reserved.</p>
      </footer>
    </div>
  );
};

export default HomePage;
