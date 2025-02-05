import { createBrowserRouter } from "react-router-dom";
import HomePage from "../../routes/HomePage";
import NotFoundPage from "../../routes/NotFoundPage";
import SearchPage from "../../routes/SearchPage";
import RootLayout from "../layout/RootLayout";

const router = createBrowserRouter([
  {
    path: "/",
    element: <RootLayout />, // Root layout for the entire app
    errorElement: <NotFoundPage />, // Fallback for unmatched routes
    children: [
      {
        index: true, // HomePage is the default route for '/'
        element: <HomePage />,
      },
      {
        path: "search", // SearchPage is accessible at '/search'
        element: <SearchPage />,
      },
    ],
  },
]);

export default router;
