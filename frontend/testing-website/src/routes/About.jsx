import { useState } from "react";
import AlexLo from "../assets/AlexLo.jpg";
import ColeOliva from "../assets/ColeOliva.jpg";
import JimHarle from "../assets/JimHarle.jpg";
import RichardJefferson from "../assets/RichardJefferson.jpg";
import ShameemAhmed from "../assets/ShameemAhmed.jpg";
import "../index.css";
import Contact from "./Contact";

const About = () => {
  const [showContact, setShowContact] = useState(false);
  return (
    <div className="max-w-4xl mx-auto py-8 px-4">
      <h1 className="text-3xl font-bold text-center mb-8">
        About the Autism Resources Database
      </h1>

      {/* Mission Section */}
      <section className="mb-10">
        <h2 className="text-2xl font-semibold mb-4">Our Mission</h2>
        <p className="text-lg mb-4">
          The Autism Resources Database (ARD) is dedicated to providing
          healthcare professionals, educators, researchers, and families with
          access to evidence-based resources and treatment papers related to
          autism spectrum disorder. Our mission is to bridge the gap between
          research and practice by making current, vetted information readily
          accessible to those who need it most.
        </p>
        <p className="text-lg">
          We believe that reliable information should be accessible to everyone
          involved in supporting individuals with autism, regardless of their
          location or resources. By centralizing research and treatment
          methodologies, we aim to promote better outcomes and quality of life
          for people on the autism spectrum.
        </p>
      </section>

      {/* Partnership Section */}
      <section className="mb-10">
        <h2 className="text-2xl font-semibold mb-4">
          Our Partnership with the Sendan Center
        </h2>
        <p className="text-lg mb-4">
          The Autism Resources Database was developed in collaboration with Jim
          Harle from the Sendan Center, a leading organization specializing in
          autism assessment, treatment, and support services. This partnership
          ensures that our database is built on clinical expertise and reflects
          the real-world needs of practitioners and families.
        </p>
        <p className="text-lg">
          The process of articles getting added to our database is thoroughly
          discussed and vetted, ensuring that all information shared through our
          platform meets high standards of scientific validity and practical
          applicability.
        </p>
      </section>

      {/* Database Content Section */}
      <section className="mb-10">
        <h2 className="text-2xl font-semibold mb-4">
          What You'll Find in Our Database
        </h2>
        <div className="grid md:grid-cols-2 gap-6">
          <div className="border border-gray-300 rounded-lg p-6">
            <h3 className="text-xl font-medium mb-3">Research Papers</h3>
            <p>
              Established studies on autism interventions, diagnostic
              approaches, and emerging treatments for mental health. Our
              collection spans from foundational research to the latest
              advancements in the field.
            </p>
          </div>
        </div>
      </section>

      {/* Team Section */}
      <section className="mb-10">
        <h2 className="text-2xl font-semibold mb-4">Our Team</h2>
        <p className="text-lg mb-6">
          The Autism Resources Database is maintained by a dedicated team of
          researchers, clinicians, and technology specialists who are passionate
          about improving access to autism resources.
        </p>
        <div className="grid md:grid-cols-3 gap-5">
          <div className="text-center">
            <img
              src={JimHarle}
              alt="Jim Harle"
              className="w-24 h-24 rounded-full object-cover mx-auto mb-4"
            />
            <h3 className="font-semibold">James Harle, MD</h3>
            <p className="text-gray-600">Child and Adolescent Psychiatrist</p>
          </div>
          <div className="text-center">
            <img
              src={ShameemAhmed}
              alt="Shameem Ahmed"
              className="w-24 h-24 rounded-full object-cover mx-auto mb-4"
            />
            <h3 className="font-semibold">Shameem Ahmed</h3>
            <p className="text-gray-600">Research Oversight and ASD Expert</p>
          </div>
          <div className="text-center">
            <img
              src={RichardJefferson}
              alt="Richard Jefferson"
              className="w-24 h-24 rounded-full object-cover mx-auto mb-4"
            />
            <h3 className="font-semibold">Richard Jefferson</h3>
            <p className="text-gray-600">Database Architect</p>
          </div>
          <div className="text-center">
            <img
              src={AlexLo}
              alt="Alex Lo"
              className="w-24 h-24 rounded-full object-cover mx-auto mb-4"
            />
            <h3 className="font-semibold">Alex Lo</h3>
            <p className="text-gray-600">API and Formatting Architect</p>
          </div>
          <div className="text-center">
            <div className="w-24 h-24 rounded-full bg-gray-300 mx-auto mb-4"></div>
            <h3 className="font-semibold">Logan Kalloway</h3>
            <p className="text-gray-600">LLM and MedBERT Architect</p>
          </div>
          <div className="text-center">
            <img
              src={ColeOliva}
              alt="Cole Oliva"
              className="w-24 h-24 rounded-full object-cover mx-auto mb-4"
            />
            <h3 className="font-semibold">Cole Oliva</h3>
            <p className="text-gray-600">Frontend and Endpoints Architect</p>
          </div>
        </div>
      </section>

      {/* Resource Vetting Section */}
      <section className="mb-10">
        <h2 className="text-2xl font-semibold mb-4">
          Our Resource Vetting Process
        </h2>
        <p className="text-lg mb-4">
          We are committed to maintaining the highest standards of accuracy and
          relevance in our database. Each resource undergoes an automated
          evaluation process before being included:
        </p>
        <ol className="list-decimal pl-6 space-y-2">
          <li>Initial screening by 3rd-party scientific databases (PubMed)</li>
          <li>
            Review by automated queries to make sure the paper is relevant to
            ASD.
          </li>
          <li>Double checking when categorizing papers to ensure accuracy</li>
          <li>
            Regular reassessment to ensure continued relevance and accuracy
          </li>
        </ol>
      </section>

      {/* Contact Section */}
      <section>
        <h2 className="text-2xl font-semibold mb-4">Get Involved</h2>
        <p className="text-lg mb-4">
          We welcome contributions from researchers, clinicians, and other
          professionals in the field. If you have resources to contribute or
          feedback on our platform, please reach out to us.
        </p>
        <p className="text-lg mb-4">
          For technical documentation, deployment guides, and developer
          resources, visit our{" "}
          <a
            href="https://tipsycanoe.github.io/Autism-Spectrum-Disorder-ASD-Treatment-Database/"
            target="_blank"
            rel="noopener noreferrer"
            className="text-blue-600 hover:text-blue-800 underline font-medium"
          >
            comprehensive documentation
          </a>
          .
        </p>
        <div className="text-center">
          <button
            onClick={() => setShowContact(true)}
            className="px-6 py-2 bg-white text-black border border-black rounded-lg hover:bg-gray-100"
          >
            Contact Us
          </button>
        </div>
      </section>

      <Contact isOpen={showContact} onClose={() => setShowContact(false)} />
    </div>
  );
};

export default About;
