import React from "react";
import "../index.css";

const About = () => {
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
          The Autism Resources Database was developed in collaboration with the
          Sendan Center, a leading organization specializing in autism
          assessment, treatment, and support services. This partnership ensures
          that our database is informed by clinical expertise and reflects the
          real-world needs of practitioners and families.
        </p>
        <p className="text-lg">
          The Sendan Center's team of specialists contributes to the vetting
          process of our resources, ensuring that all information shared through
          our platform meets high standards of scientific validity and practical
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
              Peer-reviewed studies on autism interventions, diagnostic
              approaches, and emerging therapies. Our collection spans from
              foundational research to the latest advancements in the field.
            </p>
          </div>
          <div className="border border-gray-300 rounded-lg p-6">
            <h3 className="text-xl font-medium mb-3">
              Treatment Methodologies
            </h3>
            <p>
              Detailed guides on evidence-based interventions including Applied
              Behavior Analysis (ABA), speech therapy, occupational therapy, and
              social skills training approaches.
            </p>
          </div>
          <div className="border border-gray-300 rounded-lg p-6">
            <h3 className="text-xl font-medium mb-3">Educational Resources</h3>
            <p>
              Materials for teachers and educational specialists working with
              students on the autism spectrum, including classroom strategies
              and accommodation recommendations.
            </p>
          </div>
          <div className="border border-gray-300 rounded-lg p-6">
            <h3 className="text-xl font-medium mb-3">Family Support</h3>
            <p>
              Practical guides and resources for parents and caregivers,
              addressing daily challenges, advocacy information, and approaches
              to supporting development at home.
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
        <div className="grid md:grid-cols-3 gap-6">
          <div className="text-center">
            <div className="w-32 h-32 rounded-full bg-gray-300 mx-auto mb-4"></div>
            <h3 className="font-semibold">John Smith</h3>
            <p className="text-gray-600">Research Director</p>
          </div>
          <div className="text-center">
            <div className="w-32 h-32 rounded-full bg-gray-300 mx-auto mb-4"></div>
            <h3 className="font-semibold">John Smith</h3>
            <p className="text-gray-600">Database Architect</p>
          </div>
          <div className="text-center">
            <div className="w-32 h-32 rounded-full bg-gray-300 mx-auto mb-4"></div>
            <h3 className="font-semibold">John Smith</h3>
            <p className="text-gray-600">Clinical Advisor</p>
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
          relevance in our database. Each resource undergoes a rigorous
          evaluation process before being included:
        </p>
        <ol className="list-decimal pl-6 space-y-2">
          <li>
            Initial screening by our research team for relevance and basic
            scientific merit
          </li>
          <li>
            Review by subject matter experts from various disciplines related to
            autism
          </li>
          <li>Clinical evaluation by practitioners from the Sendan Center</li>
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
        <div className="text-center">
          <button className="px-6 py-2 bg-white text-black border border-black rounded-lg hover:bg-gray-100">
            Contact Us
          </button>
        </div>
      </section>
    </div>
  );
};

export default About;
