import React from "react";
import { Link } from "react-router-dom";
import "./NavBar.css"; // Import the CSS file for NavBar styles

const NavBar = () => {
  return (
    <nav className="navbar">
      <div className = "navbar-title">ASD Database</div>
      <div className="navbar-list-wrapper">
        <ul className="navbar-list">
          <li className="navbar-item">
            <Link to="/">Home</Link>
          </li>
          <li className="navbar-item">
            <Link to="/search">Search</Link>
          </li>
          <li className="navbar-item">
            <Link to="/faq">FAQ</Link>
          </li>
          <li className="navbar-item">
            <Link to="/about">About</Link>
          </li>
        </ul>
      </div>
    </nav>
  );
};

export default NavBar;
