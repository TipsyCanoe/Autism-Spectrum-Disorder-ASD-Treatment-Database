// src/app/App.jsx
import { RouterProvider } from "react-router-dom";
import "../index.css";
import router from "./routing/router"; // Import the router


const App = () => {
  return <RouterProvider router={router} />; // Provide the router
};

export default App;
