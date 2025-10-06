# pages/2_Data_Extraction.py

import streamlit as st
from Bio import Entrez
import requests

# --- PAGE CONFIG ---
st.set_page_config(page_title="Data Extraction", page_icon="🔬", layout="wide")

# --- NCBI API Configuration ---
Entrez.email = "mo7eh7@gmail.com" 

#https://api.unpaywall.org/v2/10.1038/nature12373?email=mo7eh7@gmail.com

# --- FUNCTIONS ---
@st.cache_data
def get_doi_and_pdf_link(article_id):
    """Fetch an article's DOI from PubMed and then find its PDF link via Unpaywall."""
    # Step 1: Get the DOI from PubMed
    handle = Entrez.efetch(db="pubmed", id=article_id, rettype="fasta", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    
    doi = None
    try:
        # DOI is often in the ArticleIdList
        for item in records['PubmedArticle'][0]['PubmedData']['ArticleIdList']:
            if item.attributes['IdType'] == 'doi':
                doi = str(item)
                print(doi)
                break
    except KeyError:
        pass # No DOI found

    if not doi:
        return "DOI not found", None

    # Step 2: Query the Unpaywall API with the DOI
    unpaywall_url = f"https://api.unpaywall.org/v2/{doi}?email=mo7eh7@gmail.com"
    try:
        response = requests.get(unpaywall_url)
        response.raise_for_status()
        data = response.json()
        pdf_url = data.get('best_oa_location', {}).get('url_for_pdf')
        
        if pdf_url:
            return doi, pdf_url
        else:
            return doi, "No open-access PDF found"
            
    except requests.exceptions.RequestException:
        return doi, "Could not contact Unpaywall"

# --- USER INTERFACE ---
st.title("📝 Manara AI - Data Extraction Workspace")

# Get the list of included articles from the session state
included_articles = [
    article for article in st.session_state.get('search_results', []) 
    if article['status'] == 'included'
]

if not included_articles:
    st.warning("No articles have been marked as 'Included' in the Screening phase yet.")
    st.page_link("1_Screening.py", label="Back to Screening", icon="⬅️")
else:
    st.info(f"You have moved {len(included_articles)} articles to the data extraction phase.")
    
    # Create a dropdown to select which article to work on
    selected_title = st.selectbox(
        "Select an article to extract data from:",
        options=[a['title'] for a in included_articles]
    )
    
    # Find the full article object based on the selected title
    selected_article = next((a for a in included_articles if a['title'] == selected_title), None)

    if selected_article:
        article_id = selected_article['id']
        
        # --- Workspace Layout ---
        st.markdown("---")
        
        # Fetch DOI and PDF link
        with st.spinner("Searching for article DOI and open-access PDF..."):
            doi, pdf_link = get_doi_and_pdf_link(article_id)

        st.subheader(selected_title)
        st.caption(f"PubMed ID: {article_id} | DOI: {doi}")

        col1, col2 = st.columns(2)

        with col1:
            st.header("Full-Text PDF")
            if pdf_link and "http" in pdf_link:
                st.success("Open-access PDF found!")
                st.link_button("Open PDF in New Tab ↗️", pdf_link)
                # Display the PDF in an iframe for a more integrated feel
                st.components.v1.iframe(pdf_link, height=800, scrolling=True)
            else:
                st.error(f"PDF Status: {pdf_link}")

        with col2:
            st.header("Data Extraction Form")
            st.write("*(This is a placeholder for the data extraction template)*")
            
            # Example of what the form will look like
            st.text_input("Sample Size (n)")
            st.text_area("Patient Demographics")
            st.text_area("Primary Outcome(s)")
            st.selectbox("Risk of Bias Assessment", ["Low", "Some Concerns", "High"])
            st.button("Save Data", type="primary")