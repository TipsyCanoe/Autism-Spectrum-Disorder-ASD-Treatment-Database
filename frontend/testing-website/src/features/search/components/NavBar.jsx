import { Link } from "react-router-dom";

// top of the page navigation!!!
const NavBar = () => {
  return (
    <nav
      className="fixed top-0 w-full bg-navbar-blue shadow-md z-50 flex items-center px-5 rounded-bl-md rounded-br-md"
      role="navigation"
    >
      <div className="absolute left-5 py-2.5 text-white select-none m-0">
        {" "}
        ARD
      </div>
      <div className="flex-1 flex justify-center">
        {" "}
        <ul className="flex list-none p-0 m-0">
          {" "}
          <li className="m-0">
            {" "}
            <Link
              to="/"
              className="text-white px-5 py-3.75 block hover:bg-link-hover-blue"
            >
              Home
            </Link>
          </li>
          <li className="m-0">
            {" "}
            <Link
              to="/search"
              className="text-white px-5 py-3.75 block hover:bg-link-hover-blue"
            >
              Search
            </Link>
          </li>
          <li className="m-0">
            {" "}
            <Link
              to="/faq"
              className="text-white px-5 py-3.75 block hover:bg-link-hover-blue"
            >
              FAQ
            </Link>
          </li>
          <li className="m-0">
            {" "}
            <Link
              to="/about"
              className="text-white px-5 py-3.75 block hover:bg-link-hover-blue"
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
