import { useState } from "react";
import { Link } from "react-router-dom";
import AutismSociety from "../assets/AutismSociety.png";
import AutismSpeaksImg from "../assets/AutismSpeaksLink.jpg";
import NationalAutismAssociation from "../assets/NationalAutismAssociation.png";
import "../index.css";

const HomePage = () => {
  const [isJobRunning, setIsJobRunning] = useState(false);
  const [jobMessage, setJobMessage] = useState("");

  const handleRunJob = async () => {
    setIsJobRunning(true);
    setJobMessage("Updating database...");
    try {
      const response = await fetch("/jobapi/run-job", {
        method: "POST",
      });
      const data = await response.json();
      if (response.ok) {
        setJobMessage(
          data.message + (data.details ? ` Details: ${data.details}` : "")
        );
      } else {
        setJobMessage(
          `Error: ${data.message}` +
            (data.error ? ` Details: ${data.error}` : "")
        );
      }
    } catch (error) {
      console.error("Failed to trigger API job:", error);
      setJobMessage("Failed to trigger job. Check console for details.");
    }
    setIsJobRunning(false);
  };

  return (
    <main>
      <section className="bg-gray-100 py-12">
        <div className="container mx-auto text-center">
          <h1 className="text-4xl font-bold text-gray-800">
            Welcome to the Autism Resources Database
          </h1>
          <p className="text-gray-600 mt-4">
            A comprehensive database of resources for professionals and
            families. We are working with the{" "}
            <a
              href="https://sendancenter.com/"
              target="_blank"
              rel="noopener noreferrer"
            >
              {" "}
              <b>Sendan Center</b>
            </a>{" "}
            to provide the best resources for autism support.
          </p>
        </div>
      </section>

      {/* Search Section */}
      <section className="bg-white py-12">
        <div className="container mx-auto text-center">
          <h2 className="text-2xl font-bold text-gray-800 mb-4">
            Find the Resources You Need
          </h2>
          <div className="flex justify-center items-center space-x-4">
            <Link to="/search">
              <button className="px-6 py-2 bg-navbar-blue text-white rounded-lg hover:bg-link-hover-blue">
                Search
              </button>
            </Link>
            <button
              onClick={handleRunJob}
              disabled={isJobRunning}
              className="px-6 py-2 bg-navbar-blue text-white rounded-lg hover:bg-link-hover-blue disabled:bg-gray-400"
            >
              {isJobRunning ? "Updating..." : "Update Database"}
            </button>
          </div>
          {jobMessage && (
            <p className="mt-4 text-sm text-gray-600">{jobMessage}</p>
          )}
        </div>
      </section>

      {/* Features Section */}
      <section className="bg-white py-12">
        <div className="container mx-auto text-center">
          <h2 className="text-2xl font-bold text-gray-800 text-center mb-8">
            Features
          </h2>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h3 className="text-xl font-semibold">Extensive Database</h3>
              <p className="text-gray-600 mt-2">
                Access a wide range of resources tailored for autism support.
              </p>
            </div>
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h3 className="text-xl font-semibold">Professional Guidance</h3>
              <p className="text-gray-600 mt-2">
                Find expert advice and best practices for autism care.
              </p>
            </div>
            <div className="p-6 bg-gray-100 rounded-lg shadow">
              <h3 className="text-xl font-semibold">Family Support</h3>
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
          <h2 className="text-2xl font-bold text-gray-800 mb-8">
            What Our Users Say
          </h2>
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
          <h2 className="text-2xl font-bold text-gray-800 mb-8">
            Other Autism Resources
          </h2>
          <div className="horizontal-resources justify-center">
            <div className="flex flex-col items-center p-6 bg-gray-100 rounded-lg shadow w-80">
              <img
                src={AutismSpeaksImg}
                alt="Autism Speaks"
                className="w-16 h-16 object-cover rounded-lg mb-4"
              />

              <h3 className="text-xl font-semibold">Autism Speaks</h3>
              <p className="text-gray-600 mt-2">
                Autism Speaks is dedicated to promoting solutions for the needs
                of individuals with autism and their families.
              </p>
              <a
                href="https://www.autismspeaks.org"
                target="_blank"
                rel="noopener noreferrer"
                className="mt-4 inline-block px-6 py-2 bg-navbar-blue text-white rounded-lg hover:bg-link-hover-blue"
              >
                Visit Autism Speaks
              </a>
            </div>

            <div className="flex flex-col items-center p-6 bg-gray-100 rounded-lg shadow w-80">
              <img
                src={NationalAutismAssociation}
                alt="National Autism Association"
                className="w-16 h-16 object-cover rounded-lg mb-4"
              />

              <h3 className="text-xl font-semibold">
                National Autism Association
              </h3>
              <p className="text-gray-600 mt-2">
                The National Autism Association provides resources, education,
                advocacy, and support for families affected by autism.
              </p>
              <a
                href="https://nationalautismassociation.org"
                target="_blank"
                rel="noopener noreferrer"
                className="mt-4 inline-block px-6 py-2 bg-navbar-blue text-white rounded-lg hover:bg-link-hover-blue"
              >
                Visit National Autism Association
              </a>
            </div>

            <div className="flex flex-col items-center p-6 bg-gray-100 rounded-lg shadow w-80">
              <img
                src={AutismSociety}
                alt="Autism Society"
                className="w-16 h-16 object-cover rounded-lg mb-4"
              />

              <h3 className="text-xl font-semibold">Autism Society</h3>
              <p className="text-gray-600 mt-2">
                The Autism Society aims to improve the lives of all affected by
                autism through advocacy, education, and support.
              </p>
              <a
                href="https://www.autism-society.org"
                target="_blank"
                rel="noopener noreferrer"
                className="mt-4 inline-block px-6 py-2 bg-navbar-blue text-white rounded-lg hover:bg-link-hover-blue"
              >
                Visit Autism Society
              </a>
            </div>
          </div>
        </div>
      </section>

      <footer className="bg-gray-800 text-white text-center py-4">
        <p>&copy; 2025 Autism Resources Database. All rights reserved.</p>
      </footer>
    </main>
  );
};

export default HomePage;
