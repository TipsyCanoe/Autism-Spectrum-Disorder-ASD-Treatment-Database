import { useState } from "react";
import { Link } from "react-router-dom";
import AutismSociety from "../assets/AutismSociety.png";
import NationalAutismAssociation from "../assets/NationalAutismAssociation.png";
import "../index.css";

const HomePage = () => {
  const [isJobRunning, setIsJobRunning] = useState(false);
  const [jobMessage, setJobMessage] = useState("");

  const handleRunJob = async () => {
    setIsJobRunning(true);
    setJobMessage("Updating database...");
    try {
      const apiUrl = process.env.REACT_APP_NODE_API_URL || "http://localhost:5001";
      const response = await fetch(`${apiUrl}/api/run-job`, {
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
      <section className="bg-gray-100 py-16">
        <div className="container mx-auto text-center px-4">
          <h1 className="text-4xl lg:text-5xl font-bold text-gray-800">
            Welcome to STAR
          </h1>
          <p className="card-heading text-gray-700 mt-3">
            Sendan Tools for Autism Resources
          </p>
          <p className="body-text text-gray-600 mt-6 max-w-3xl mx-auto leading-relaxed">
            A comprehensive database of peer-reviewed research for professionals and
            families. We are working with the{" "}
            <a
              href="https://sendancenter.com/"
              target="_blank"
              rel="noopener noreferrer"
              className="text-navbar-blue hover:underline font-semibold"
            >
              Sendan Center
            </a>{" "}
            to aggregate and organize autism treatment studies from PubMed.
          </p>
        </div>
      </section>

      {/* Search Section */}
      <section className="bg-white py-16">
        <div className="container mx-auto text-center px-4">
          <h2 className="section-heading text-gray-800 text-center">
            Find the Research You Need
          </h2>
          <div className="flex flex-col items-center space-y-4">
            <div className="flex flex-wrap justify-center items-center gap-4">
              <Link to="/search">
                <button className="btn-primary">
                  Search Database
                </button>
              </Link>
              <button
                onClick={handleRunJob}
                disabled={isJobRunning}
                className="btn-primary disabled:bg-gray-400"
                title="Fetches the latest studies from PubMed"
              >
                {isJobRunning ? "Updating..." : "Update Database"}
              </button>
            </div>
            <p className="body-text-sm text-gray-500 max-w-2xl">
              The database automatically pulls the latest autism treatment studies from PubMed. 
              Click "Update Database" to manually trigger a refresh.
            </p>
          </div>
          {jobMessage && (
            <p className="mt-4 body-text-sm text-gray-600">{jobMessage}</p>
          )}
        </div>
      </section>

      {/* Statistics Section */}
      <section className="bg-white py-16 border-y border-gray-200">
        <div className="container mx-auto px-4">
          <h2 className="section-heading text-gray-800 text-center">
            Database Overview
          </h2>
          <div className="grid grid-cols-1 md:grid-cols-4 gap-8 max-w-6xl mx-auto">
            <div className="text-center">
              <div className="text-4xl lg:text-5xl font-bold text-navbar-blue mb-3">1000+</div>
              <p className="body-text text-gray-700 font-medium">Research Studies</p>
              <p className="body-text-sm text-gray-600 mt-2">Peer-reviewed publications</p>
            </div>
            <div className="text-center">
              <div className="text-4xl lg:text-5xl font-bold text-navbar-blue mb-3">150+</div>
              <p className="body-text text-gray-700 font-medium">Treatment Types</p>
              <p className="body-text-sm text-gray-600 mt-2">Across multiple categories</p>
            </div>
            <div className="text-center">
              <div className="text-4xl lg:text-5xl font-bold text-navbar-blue mb-3">25+</div>
              <p className="body-text text-gray-700 font-medium">Years Covered</p>
              <p className="body-text-sm text-gray-600 mt-2">Recent and historical data</p>
            </div>
            <div className="text-center">
              <div className="text-4xl lg:text-5xl font-bold text-navbar-blue mb-3">10+</div>
              <p className="body-text text-gray-700 font-medium">Search Filters</p>
              <p className="body-text-sm text-gray-600 mt-2">Refined search options</p>
            </div>
          </div>
        </div>
      </section>

      {/* Features Section */}
      <section className="bg-gray-50 py-16">
        <div className="container mx-auto text-center px-4">
          <h2 className="section-heading text-gray-800 text-center">
            What We Offer
          </h2>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8 lg:gap-10">
            <Link to="/search" className="block">
              <div className="p-8 bg-gray-100 rounded-lg shadow hover:shadow-lg hover:bg-gray-200 transition-all cursor-pointer h-full">
                <h3 className="card-heading text-navbar-blue">Extensive Database</h3>
                <p className="body-text-sm text-gray-700 mt-4 leading-relaxed">
                  Access a wide range of peer-reviewed studies and treatment research tailored for autism support.
                </p>
                <span className="body-text-sm text-navbar-blue mt-6 inline-block font-medium">
                  Explore Database →
                </span>
              </div>
            </Link>
            <Link to="/about" className="block">
              <div className="p-8 bg-gray-100 rounded-lg shadow hover:shadow-lg hover:bg-gray-200 transition-all cursor-pointer h-full">
                <h3 className="card-heading text-navbar-blue">Curated from PubMed</h3>
                <p className="body-text-sm text-gray-700 mt-4 leading-relaxed">
                  Studies sourced from PubMed's scientific database and organized in partnership with the Sendan Center.
                </p>
                <span className="body-text-sm text-navbar-blue mt-6 inline-block font-medium">
                  Learn More →
                </span>
              </div>
            </Link>
            <Link to="/search" className="block">
              <div className="p-8 bg-gray-100 rounded-lg shadow hover:shadow-lg hover:bg-gray-200 transition-all cursor-pointer h-full">
                <h3 className="card-heading text-navbar-blue">Evidence-Based Resources</h3>
                <p className="body-text-sm text-gray-700 mt-4 leading-relaxed">
                  Find validated treatment approaches and intervention studies to support informed decisions.
                </p>
                <span className="body-text-sm text-navbar-blue mt-6 inline-block font-medium">
                  Start Searching →
                </span>
              </div>
            </Link>
          </div>
        </div>
      </section>

      {/* How It Works Section */}
      <section className="bg-white py-16">
        <div className="container mx-auto px-4">
          <h2 className="section-heading text-gray-800 text-center">
            How It Works
          </h2>
          <div className="max-w-4xl mx-auto space-y-8">
            <div className="flex items-start gap-6">
              <div className="flex-shrink-0 w-16 h-16 bg-navbar-blue text-white rounded-full flex items-center justify-center text-2xl font-bold">
                1
              </div>
              <div>
                <h3 className="card-heading text-gray-800 mb-3">Search or Browse</h3>
                <p className="body-text-sm text-gray-700 leading-relaxed">
                  Use our advanced filters to search by treatment type, age group, outcome measures, or publication year. You can also browse all studies without filters.
                </p>
              </div>
            </div>
            
            <div className="flex items-start gap-6">
              <div className="flex-shrink-0 w-16 h-16 bg-navbar-blue text-white rounded-full flex items-center justify-center text-2xl font-bold">
                2
              </div>
              <div>
                <h3 className="card-heading text-gray-800 mb-3">Review Results</h3>
                <p className="body-text-sm text-gray-700 leading-relaxed">
                  Explore study summaries including treatment details, participant information, outcomes, and statistical measures. All studies are sourced from PubMed.
                </p>
              </div>
            </div>
            
            <div className="flex items-start gap-6">
              <div className="flex-shrink-0 w-16 h-16 bg-navbar-blue text-white rounded-full flex items-center justify-center text-2xl font-bold">
                3
              </div>
              <div>
                <h3 className="card-heading text-gray-800 mb-3">Make Informed Decisions</h3>
                <p className="body-text-sm text-gray-700 leading-relaxed">
                  Access direct links to full publications and use the evidence-based information to support your research or decision-making process.
                </p>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Other Resources Section */}
      <section className="bg-gray-50 py-16">
        <div className="container mx-auto text-center px-4">
          <h2 className="section-heading text-gray-800 text-center">
            Other Autism Resources
          </h2>
          <div className="horizontal-resources justify-center">
            <div className="flex flex-col items-center p-8 bg-gray-100 rounded-lg shadow w-80 lg:w-96">
              <img
                src={NationalAutismAssociation}
                alt="National Autism Association"
                className="w-20 h-20 object-cover rounded-lg mb-6"
              />

              <h3 className="card-heading mb-4">
                National Autism Association
              </h3>
              <p className="body-text-sm text-gray-700 mt-2 leading-relaxed">
                The National Autism Association provides resources, education,
                advocacy, and support for families affected by autism.
              </p>
              <a
                href="https://nationalautismassociation.org"
                target="_blank"
                rel="noopener noreferrer"
                className="btn-primary mt-6"
              >
                Visit National Autism Association
              </a>
            </div>

            <div className="flex flex-col items-center p-8 bg-gray-100 rounded-lg shadow w-80 lg:w-96">
              <img
                src={AutismSociety}
                alt="Autism Society"
                className="w-20 h-20 object-cover rounded-lg mb-6"
              />

              <h3 className="card-heading mb-4">Autism Society</h3>
              <p className="body-text-sm text-gray-700 mt-2 leading-relaxed">
                The Autism Society aims to improve the lives of all affected by
                autism through advocacy, education, and support.
              </p>
              <a
                href="https://www.autism-society.org"
                target="_blank"
                rel="noopener noreferrer"
                className="btn-primary mt-6"
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
