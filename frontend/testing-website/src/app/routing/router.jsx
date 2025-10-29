import { createBrowserRouter } from "react-router-dom";
import About from "../../routes/About";
import FAQ from "../../routes/FAQ";
import HomePage from "../../routes/HomePage";
import NotFoundPage from "../../routes/NotFoundPage";
import SearchPage from "../../routes/SearchPage";
import RootLayout from "../layout/RootLayout";

const router = createBrowserRouter([
  {
    path: "/",
    element: <RootLayout />,
    errorElement: <NotFoundPage />,
    children: [
      { index: true, element: <HomePage /> },
      {
        path: "search",
        element: <SearchPage />,
      },
      { path: "FAQ", element: <FAQ /> },
      { path: "about", element: <About /> },
    ],
  },
]);

export default router;
