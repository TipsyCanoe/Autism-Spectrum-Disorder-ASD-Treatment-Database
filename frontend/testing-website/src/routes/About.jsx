import { useState } from "react";
import AlexLo from "../assets/AlexLo.jpg";
import ColeOliva from "../assets/ColeOliva.jpg";
import JimHarle from "../assets/JimHarle.jpg";
import LoganKalloway from "../assets/LoganKalloway.jpg";
import RichardJefferson from "../assets/RichardJefferson.jpg";
import ShameemAhmed from "../assets/ShameemAhmed.jpg";
import "../index.css";
import Contact from "./Contact";

const About = () => {
  const [showContact, setShowContact] = useState(false);
  return (
    <div className="max-w-4xl mx-auto py-8 px-4">
      <h1 className="page-heading">
        About the Autism Resources Database
      </h1>

      {/* Mission Section */}
      <section className="mb-10">
        <h2 className="section-heading">Our Mission</h2>
        <p className="body-text mb-4">
          The Autism Resources Database (ARD) is dedicated to providing
          healthcare professionals, educators, researchers, and families with
          access to evidence-based resources and treatment papers related to
          autism spectrum disorder. Our mission is to bridge the gap between
          research and practice by making current, vetted information readily
          accessible to those who need it most.
        </p>
        <p className="body-text">
          We believe that reliable information should be accessible to everyone
          involved in supporting individuals with autism, regardless of their
          location or resources. By centralizing research and treatment
          methodologies, we aim to promote better outcomes and quality of life
          for people on the autism spectrum.
        </p>
      </section>

      {/* Partnership Section */}
      <section className="mb-10">
        <h2 className="section-heading">
          Our Partnership with the Sendan Center
        </h2>
        <p className="body-text mb-4">
          The Autism Resources Database was developed in collaboration with Jim
          Harle from the Sendan Center, a leading organization specializing in
          autism assessment, treatment, and support services. This partnership
          ensures that our database is built on clinical expertise and reflects
          the real-world needs of practitioners and families.
        </p>
        <p className="body-text">
          The process of articles getting added to our database is thoroughly
          discussed and vetted, ensuring that all information shared through our
          platform meets high standards of scientific validity and practical
          applicability.
        </p>
      </section>

      {/* Database Content Section */}
      <section className="mb-10">
        <h2 className="section-heading">
          What You'll Find in Our Database
        </h2>
        <div className="grid md:grid-cols-2 gap-6">
          <div className="border border-gray-300 rounded-lg p-6">
            <h3 className="card-heading mb-3">Research Papers</h3>
            <p className="body-text">
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
        <h2 className="section-heading">Our Team</h2>
        <p className="body-text mb-6">
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
            <img
              src={LoganKalloway}
              alt="Logan Kalloway"
              className="w-24 h-24 rounded-full object-cover mx-auto mb-4"
            />
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
        <h2 className="section-heading">
          Our Resource Vetting Process
        </h2>
        <p className="body-text mb-4">
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
        <h2 className="section-heading">Get Involved</h2>
        <p className="body-text mb-4">
          We welcome contributions from researchers, clinicians, and other
          professionals in the field. If you have resources to contribute or
          feedback on our platform, please reach out to us.
        </p>
        <p className="body-text mb-4">
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
            className="btn-standard bg-white text-black border border-black hover:bg-gray-100"
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
