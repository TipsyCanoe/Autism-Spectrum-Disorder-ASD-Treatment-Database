# Autism-Spectrum-Disorder-ASD-Treatment-Database

Vision: To enhance the mental health of individuals with Autism Spectrum Disorder (ASD)
and their families by synthesizing psychiatric treatment knowledge for healthcare
professionals, patients, and families

The goal of this project is to create a database and front-end representation for families and
healthcare professionals looking to know more about autism. Resources will be pulled from various
published studies and articles.

If you want to run the frontend-barebones testing website locally, get into the frontend directory, and then the testing-website directory inside. 
Then, run ```npm install```. After that, running ```npm start``` should be sufficient to deploy locally. if you ctrl-z to 
end the process in the terminal, the process still may be running on port 3000 (I think it's specified). run ```fuser -k 3000/tcp``` 
which should kill the process. 

## PubMed API Extraction
Source: https://youtu.be/sGC66q45BX4

The code was cloned from a public Github repository, obtained from the above Youtube video.
No email needed for now, as I doubt we are using more than 3 requests per second.

Dependencies (if using Python 3, just add 3 after pip or python): 
    ```sh
    pip install pandas biopython
    ```
Run the python program, and output is in PubMed_results.xlsx.
For full README, go to the github listed in the video.
Full-text download is possible, in xml, txt, and pdf format I think, but I haven't gotten that yet.
