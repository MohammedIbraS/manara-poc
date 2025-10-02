# main.py

import streamlit as st
import datetime
from Bio import Entrez
from sentence_transformers import SentenceTransformer, util # NEW: Import sentence-transformers

# --- PAGE CONFIG ---
st.set_page_config(page_title="Manara Screening", page_icon="🔬", layout="wide")

# --- NCBI API Configuration ---
Entrez.email = "your.email@example.com" 

# --- MODEL LOADING ---
@st.cache_resource
def load_relevance_model():
    """Load the Sentence Transformer model for calculating relevance."""
    return SentenceTransformer('all-MiniLM-L6-v2')

model = load_relevance_model()

# --- INITIALIZE SESSION STATE ---
if 'search_results' not in st.session_state:
    st.session_state.search_results = []
if 'error' not in st.session_state:
    st.session_state.error = None

# --- FUNCTIONS ---
def construct_pico_query(p=None, i=None, c=None, o=None, saudi_filter=False, year_range=None):
    """Construct a high-precision PubMed query from all components."""
    terms = []
    if p: terms.append(f'"{p}"[Title/Abstract]')
    if i: terms.append(f'"{i}"[Title/Abstract]')
    if c: terms.append(f'"{c}"[Title/Abstract]')
    if o: terms.append(f'"{o}"[Title/Abstract]')

    if saudi_filter:
        terms.append(f'("Saudi Arabia"[MeSH Terms] OR "Saudi Arabia"[Title/Abstract])')
    
    if year_range:
        start_year, end_year = year_range
        terms.append(f'("{start_year}"[Date - Publication] : "{end_year}"[Date - Publication])')

    if not terms:
        return None
    return " AND ".join(terms)

def get_abstract(article_id):
    """Fetch the abstract for a single PubMed article ID using Biopython."""
    try:
        handle = Entrez.efetch(db="pubmed", id=article_id, rettype="abstract", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        abstract = record['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
        return str(abstract)
    except (IndexError, KeyError):
        return None # Return None if abstract is not found
    except Exception:
        return None

def search_pubmed(query, pico_query_text, max_results=10):
    """Search PubMed, calculate relevance, and update the session state."""
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=str(max_results), sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        
        if not id_list:
            st.session_state.search_results = []
            st.session_state.error = "No articles found with that specific combination of filters."
            return

        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="fasta", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        results = []
        abstracts_to_score = []
        for record in records['PubmedArticle']:
            article_id = record['MedlineCitation']['PMID']
            title = record['MedlineCitation']['Article']['ArticleTitle']
            abstract = get_abstract(article_id)
            
            if abstract: # Only process articles that have an abstract
                results.append({"title": title, "id": article_id, "status": "pending", "abstract": abstract})
                abstracts_to_score.append(abstract)

        # --- NEW: AI Relevance Scoring ---
        if abstracts_to_score:
            query_embedding = model.encode(pico_query_text, convert_to_tensor=True)
            abstract_embeddings = model.encode(abstracts_to_score, convert_to_tensor=True)
            cosine_scores = util.cos_sim(query_embedding, abstract_embeddings)
            
            # Add scores to the results
            for i, article in enumerate(results):
                article['relevance_score'] = cosine_scores[0][i].item() # Get score as a float
            
            # Sort the results by the new relevance score, highest first
            results.sort(key=lambda x: x['relevance_score'], reverse=True)

        st.session_state.search_results = results
        st.session_state.error = None

    except Exception as e:
        st.session_state.search_results = []
        st.session_state.error = f"An API error occurred: {e}"


def set_status(article_id, new_status):
    """Find the article in session_state and update its status."""
    for article in st.session_state.search_results:
        if article['id'] == article_id:
            article['status'] = new_status
            break

# --- USER INTERFACE ---
st.title("🔬 Manara AI - Screening Workspace")

with st.sidebar:
    st.header("Search Controls")
    p = st.text_input("**P** - Population/Problem")
    i = st.text_input("**I** - Intervention")
    c = st.text_input("**C** - Comparison (Optional)")
    o = st.text_input("**O** - Outcome (Optional)")

    current_year = datetime.date.today().year
    year_range = st.slider("Publication Year Range", min_value=1950, max_value=current_year, value=(current_year - 10, current_year))
    
    saudi_filter = st.checkbox("Filter for Saudi Arabia only", value=True)

    if st.button("Search PubMed", type="primary"):
        pico_components = [p, i, c, o]
        pico_query_text = " ".join(filter(None, pico_components)) # Create a single string from PICO inputs
        
        query = construct_pico_query(p, i, c, o, saudi_filter=saudi_filter, year_range=year_range)
        if query:
            with st.spinner("Searching PubMed & calculating relevance scores..."):
                search_pubmed(query, pico_query_text)
        else:
            st.warning("Please enter at least one PICO term.")
    
    st.divider()
    st.header("Screening Progress")
    if st.session_state.search_results:
        included_count = sum(1 for a in st.session_state.search_results if a['status'] == 'included')
        excluded_count = sum(1 for a in st.session_state.search_results if a['status'] == 'excluded')
        pending_count = sum(1 for a in st.session_state.search_results if a['status'] == 'pending')
        
        st.metric("✅ Included", included_count)
        st.metric("❌ Excluded", excluded_count)
        st.metric("🤔 Pending", pending_count)

st.header("Search Results")

if st.session_state.error:
    st.error(st.session_state.error)
elif not st.session_state.search_results:
    st.info("Enter your search criteria in the sidebar and click 'Search PubMed' to begin.")
else:
    for article in st.session_state.search_results:
        # --- NEW: Display the relevance score with the title ---
        score_percentage = article.get('relevance_score', 0) * 100
        display_title = f"**{score_percentage:.1f}% Match** | {article['title']}"
        
        with st.expander(f"*{article['status'].upper()}*: {display_title}"):
            st.markdown(article.get('abstract', 'Abstract not available.'))
            st.divider()
            
            col1, col2 = st.columns(2)
            with col1:
                st.button("✅ Include", key=f"include_{article['id']}", on_click=set_status, args=(article['id'], '✅ included'), use_container_width=True)
            with col2:
                st.button("❌ Exclude", key=f"exclude_{article['id']}", on_click=set_status, args=(article['id'], '❌ excluded'), use_container_width=True)