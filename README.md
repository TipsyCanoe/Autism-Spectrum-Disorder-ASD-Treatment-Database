# Autism-Spectrum-Disorder-ASD-Treatment-Database

Vision: To enhance the mental health of individuals with Autism Spectrum Disorder (ASD)
and their families by synthesizing psychiatric treatment knowledge for healthcare
professionals, patients, and families

The goal of this project is to create a database and front-end representation for families and
healthcare professionals looking to know more about autism. Resources will be pulled from various
published studies and articles.

If you want to run the frontend testing website locally, get into the frontend directory, and then the testing-website directory inside (``cd frontend/testing-website``)
Then, run ``npm install``. After that, running ``npm start`` should be sufficient to deploy locally. if you ctrl-z to
end the process in the terminal, the process still may be running on port 3000 (I think it's specified). run ``fuser -k 3000/tcp`` or ``npx kill-port 3000``
which should kill the process.

To run the dummy backend, temporarily, in the backend directory, run ``npm run dev``. In a different terminal, for the frontend, run the ``npm start`` command in the correct directory, as detailed above. The backend should run locally on port 5001.

For testing purposes, run ``npm test`` inside the frontend/testing-website/ directory. To visualize code coverage with the tests, run ``npm test -- --coverage --watchAll=false``

## PubMed API Extraction

Source: https://youtu.be/sGC66q45BX4

I started the code base from a public Github repository, referenced from the above Youtube video.
To increase the number of requests I can make per second to the API, I have included my email. I have not included my API key, as Github does not like exposed secrets. If you want to make faster requests, enter your own email and API key, freely obtainable through signing in with your Western account.

I have two extraction files currently; one based on the query Jim sent, and one I crafted on my own with the database Jim sent as well. The query based one is ``pubmed_API_data.py``, while the database one is `pubmed_API_ASD_data.py`. Running each will output an XML sheet of `pubmed_papers_info.xlsx` and `Adjusted_ASD_Sheet_Info_V2.xlsx` respectively.

`convert_excel_to_json.py` is as the name implies; it converts the excel spreadsheet to a JSON. Currently, it expects the two excel sheets listed above.

Planned additions:
- Automatic querying based on keywords using files (currently working + unfinished)
- Obtain full-text; if not for every, then at least some (have an idea, unimplemented)
- Different queries to meet all the trial types and papers for customer request (currently working on with customer)

Dependencies (if using Python 3, just add 3 after pip or python): ``sh pip install pandas biopython``


### Using the medBERT/LLM
1. Create a virtual environment if it's your first time running
python3 -m venv venv
2. Activate the environment
source venv/bin/activate
3. Install the packages
pip install -r requirements.txt
4. You should now be able to run the LLM and MedBERT

### Using the backend
1. Create a virtual environment if it's your first time running
python3 -m venv venv
2. Activate the environment
source venv/bin/activate
3. Install the packages (in backend if you don't want it to take a while)
pip install -r requirements.txt
4. run app.py and the frontend on seperate terminals