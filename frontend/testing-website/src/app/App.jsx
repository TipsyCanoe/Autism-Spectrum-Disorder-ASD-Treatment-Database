// src/app/App.jsx
import { RouterProvider } from "react-router-dom";
import "../index.css";
import router from "./routing/router";

const App = () => {
  return <RouterProvider router={router} />;
};

export default App;
