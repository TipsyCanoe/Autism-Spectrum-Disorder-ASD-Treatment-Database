import React from "react";
import { Link } from "react-router-dom";
import AutismSpeaksImg from "../assets/AutismSpeaksLink.jpg";
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
        </div>
      </section>

      {/* Search Section */}
      <section className="bg-white py-16">
        <div className="container mx-auto text-center">
          <h3 className="text-2xl font-bold text-gray-800 mb-4">
            Find the Resources You Need
          </h3>
          <Link to="/search">
            <button className="mt-4 px-6 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700">
              Search
            </button>
          </Link>
        </div>
      </section>

      {/* Features Section */}
      <section className="bg-white py-16">
        <div className="container mx-auto text-center">
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

      {/* Testimonials Section */}
      <section className="bg-gray-100 py-16">
        <div className="container mx-auto text-center">
          <h3 className="text-2xl font-bold text-gray-800 mb-8">
            What Our Users Say
          </h3>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
            <div className="p-6 bg-white rounded-lg shadow">
              <p className="text-gray-600">
                "This database has been a lifesaver for our family. We found so
                many helpful resources!"
              </p>
              <p className="text-gray-800 mt-4 font-semibold">- Jane Doe</p>
            </div>
            <div className="p-6 bg-white rounded-lg shadow">
              <p className="text-gray-600">
                "As a professional, I rely on this database for the latest
                research and best practices."
              </p>
              <p className="text-gray-800 mt-4 font-semibold">- John Smith</p>
            </div>
          </div>
        </div>
      </section>

      {/* Other Resources Section */}
      <section className="bg-white py-16">
        <div className="container mx-auto text-center">
          <h3 className="text-2xl font-bold text-gray-800 mb-8">
            Other Autism Resources
          </h3>
          <div className="flex flex-wrap justify-center gap-8">
            <div className="flex flex-col items-center p-6 bg-gray-100 rounded-lg shadow w-80">
              <img
                src={AutismSpeaksImg}
                alt="Autism Speaks"
                className="w-24 h-24 object-cover rounded-lg mb-4"
              />
              <h4 className="text-xl font-semibold">Autism Speaks</h4>
              <p className="text-gray-600 mt-2">
                Autism Speaks is dedicated to promoting solutions for the needs
                of individuals with autism and their families.
              </p>
              <a
                href="https://www.autismspeaks.org"
                target="_blank"
                rel="noopener noreferrer"
                className="mt-4 inline-block px-6 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700"
              >
                Visit Autism Speaks
              </a>
            </div>
            <div className="flex flex-col items-center p-6 bg-gray-100 rounded-lg shadow w-80">
              <img
                src="https://nationalautismassociation.org/wp-content/uploads/2019/02/NAA-Logo.png"
                alt="National Autism Association"
                className="w-24 h-24 object-cover rounded-lg mb-4"
              />
              <h4 className="text-xl font-semibold">
                National Autism Association
              </h4>
              <p className="text-gray-600 mt-2">
                The National Autism Association provides resources, education,
                advocacy, and support for families affected by autism.
              </p>
              <a
                href="https://nationalautismassociation.org"
                target="_blank"
                rel="noopener noreferrer"
                className="mt-4 inline-block px-6 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700"
              >
                Visit National Autism Association
              </a>
            </div>
            <div className="flex flex-col items-center p-6 bg-gray-100 rounded-lg shadow w-80">
              <img
                src="https://www.autism-society.org/wp-content/uploads/2019/04/Autism-Society-Logo.png"
                alt="Autism Society"
                className="w-24 h-24 object-cover rounded-lg mb-4"
              />
              <h4 className="text-xl font-semibold">Autism Society</h4>
              <p className="text-gray-600 mt-2">
                The Autism Society aims to improve the lives of all affected by
                autism through advocacy, education, and support.
              </p>
              <a
                href="https://www.autism-society.org"
                target="_blank"
                rel="noopener noreferrer"
                className="mt-4 inline-block px-6 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700"
              >
                Visit Autism Society
              </a>
            </div>
          </div>
        </div>
      </section>

      {/* Call to Action Section */}
      <section className="bg-blue-600 py-16 text-white text-center">
        <div className="container mx-auto">
          <h3 className="text-2xl font-bold mb-4">Join Our Community</h3>
          <p className="mb-8">
            Sign up to receive updates and access exclusive resources.
          </p>
          <button className="px-6 py-2 bg-white text-blue-600 rounded-lg hover:bg-gray-200">
            Sign Up Now
          </button>
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
