import React from "react";
import { Outlet } from "react-router-dom";
import Navbar from "../../features/search/components/NavBar";
import "./RootLayout.css";

const RootLayout = () => {
  return (
    <div>
      {" "}
      {/* Adjust padding to match navbar height */}
      <Navbar />
      <div
        style={{
          paddingTop: "60px",
          paddingLeft: "20px",
          paddingRight: "20px",
        }}
      >
        <Outlet /> {/* This will render the child routes */}
      </div>
    </div>
  );
};

export default RootLayout;
