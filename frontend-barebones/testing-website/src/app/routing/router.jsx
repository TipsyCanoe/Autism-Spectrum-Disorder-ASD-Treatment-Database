import { createBrowserRouter } from "react-router-dom";
import HomePage from "../../routes/HomePage";
import NotFoundPage from "../../routes/NotFoundPage";
import SearchPage from "../../routes/SearchPage";
import RootLayout from "../layout/RootLayout";
import SearchResults from "../../routes/SearchResults";
import About from "../../routes/About";
import FAQ from "../../routes/FAQ";

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
        children: [
          {
            path: "results", // SearchResults is accessible at '/search/results'
            element: <SearchResults />,
          },
        ],
      },
      {
        path: "FAQ", // FAQPage is accessible at '/FAQ'
        element: <FAQ />,
      },
      {
        path: "about", // AboutPage is accessible at '/about'
        element: <About />,
      }
    ],
  },
]);

export default router;
