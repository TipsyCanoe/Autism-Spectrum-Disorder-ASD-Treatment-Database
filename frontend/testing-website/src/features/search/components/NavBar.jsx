import React from "react";
import { Link } from "react-router-dom";

const NavBar = () => {
  return (
    <nav
      className="fixed top-0 w-full bg-navbar-blue shadow-md z-50 flex items-center px-5 rounded-bl-md rounded-br-md" // Added Tailwind classes for .navbar
      role="navigation"
    >
      <div className="absolute left-5 py-2.5 text-white select-none m-0">
        {" "}
        {/* Added Tailwind classes for .navbar-title */}
        ARD
      </div>
      <div className="flex-1 flex justify-center">
        {" "}
        {/* Added Tailwind classes for .navbar-list-wrapper */}
        <ul className="flex list-none p-0 m-0">
          {" "}
          {/* Added Tailwind classes for .navbar-list */}
          <li className="m-0">
            {" "}
            {/* Added Tailwind classes for .navbar-item */}
            <Link
              to="/"
              className="text-white px-5 py-3.75 block hover:bg-link-hover-blue" // Added Tailwind classes for .navbar-item a and a:hover
            >
              Home
            </Link>
          </li>
          <li className="m-0">
            {" "}
            {/* Added Tailwind classes for .navbar-item */}
            <Link
              to="/search"
              className="text-white px-5 py-3.75 block hover:bg-link-hover-blue" // Added Tailwind classes for .navbar-item a and a:hover
            >
              Search
            </Link>
          </li>
          <li className="m-0">
            {" "}
            {/* Added Tailwind classes for .navbar-item */}
            <Link
              to="/faq"
              className="text-white px-5 py-3.75 block hover:bg-link-hover-blue" // Added Tailwind classes for .navbar-item a and a:hover
            >
              FAQ
            </Link>
          </li>
          <li className="m-0">
            {" "}
            {/* Added Tailwind classes for .navbar-item */}
            <Link
              to="/about"
              className="text-white px-5 py-3.75 block hover:bg-link-hover-blue" // Added Tailwind classes for .navbar-item a and a:hover
            >
              About
            </Link>
          </li>
        </ul>
      </div>
    </nav>
  );
};

export default NavBar;
