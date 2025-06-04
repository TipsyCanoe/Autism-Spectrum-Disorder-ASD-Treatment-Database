const Contact = ({ isOpen, onClose }) => {
  if (!isOpen) return null;

  // the contact us thingie
  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
      <div className="bg-white p-6 rounded-lg w-full max-w-md">
        <h2 className="text-2xl font-bold text-black mb-4">Contact Us</h2>

        <form action="https://formspree.io/f/xpwpbeag" method="POST">
          <div className="mb-4">
            <label htmlFor="contact-name" className="block text-black mb-2">
              Name
            </label>
            <input
              type="text"
              id="contact-name"
              name="name"
              className="w-full p-2 border border-gray-300 rounded"
              required
            />
          </div>
          <div className="mb-4">
            <label htmlFor="contact-email" className="block text-black mb-2">
              Email
            </label>
            <input
              type="email"
              id="contact-email"
              name="email"
              className="w-full p-2 border border-gray-300 rounded"
              required
            />
          </div>
          <div className="mb-4">
            <label htmlFor="contact-message" className="block text-black mb-2">
              Message
            </label>
            <textarea
              id="contact-message"
              name="message"
              className="w-full p-2 border border-gray-300 rounded h-32"
              required
            ></textarea>
          </div>
          <div className="flex justify-end">
            <button
              type="button"
              onClick={onClose}
              className="mr-2 px-4 py-2 bg-white text-black border border-black rounded"
            >
              Cancel
            </button>
            <button
              type="submit"
              className="px-4 py-2 bg-black text-white rounded"
            >
              Send
            </button>
          </div>
        </form>
      </div>
    </div>
  );
};

export default Contact;
