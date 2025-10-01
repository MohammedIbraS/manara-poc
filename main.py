# main.py

import streamlit as st
import requests

# --- PAGE CONFIG ---
st.set_page_config(page_title="Manara Screening", page_icon="🔬", layout="wide")

# --- INITIALIZE SESSION STATE ---
# We now need a more structured state to hold articles and their screening status
if 'search_results' not in st.session_state:
    st.session_state.search_results = [] # A list of dicts: {'title': str, 'id': str, 'status': str}
if 'error' not in st.session_state:
    st.session_state.error = None

# --- FUNCTIONS ---
def construct_pico_query(p=None, i=None, c=None, o=None, saudi_filter=False):
    """Construct a high-precision PubMed query from PICO components."""
    terms = []
    if p: terms.append(f'"{p}"[Title/Abstract]')
    if i: terms.append(f'"{i}"[Title/Abstract]')
    if c: terms.append(f'"{c}"[Title/Abstract]')
    if o: terms.append(f'"{o}"[Title/Abstract]')

    if saudi_filter:
        terms.append(f'("Saudi Arabia"[MeSH Terms] OR "Saudi Arabia"[Title/Abstract])')
    
    if not terms:
        return None
    return " AND ".join(terms)

def search_pubmed(query, max_results=10):
    """Search PubMed and update the session state with a list of articles."""
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {'db': 'pubmed', 'term': query, 'retmax': str(max_results), 'retmode': 'json', 'sort': 'relevance'}
    
    try:
        response = requests.get(esearch_url, params=params)
        response.raise_for_status()
        data = response.json()
        id_list = data['esearchresult']['idlist']
        
        if not id_list:
            st.session_state.search_results = []
            st.session_state.error = "No articles found with that specific PICO combination."
            return

        esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {'db': 'pubmed', 'id': ",".join(id_list), 'retmode': 'json'}
        response = requests.get(esummary_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        # NEW: Populate search_results with a list of dictionaries, each with a status
        results = []
        for uid in data['result']['uids']:
            results.append({
                "title": data['result'][uid]['title'],
                "id": uid,
                "status": "pending" # Initial status for all articles
            })
        st.session_state.search_results = results
        st.session_state.error = None

    except requests.exceptions.RequestException as e:
        st.session_state.search_results = []
        st.session_state.error = f"An API error occurred: {e}"
    except KeyError:
        st.session_state.search_results = []
        st.session_state.error = "Could not parse PubMed response."

def get_abstract(article_id):
    """Fetch the abstract for a single PubMed article ID."""
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {'db': 'pubmed', 'id': article_id, 'retmode': 'xml', 'rettype': 'abstract'}
    
    try:
        response = requests.get(efetch_url, params=params)
        response.raise_for_status()
        full_text = response.text
        abstract = full_text.split('<AbstractText>')[1].split('</AbstractText>')[0]
        return abstract
    except IndexError:
        return "Abstract not available in a standard format."
    except requests.exceptions.RequestException:
        return "Could not fetch abstract due to a network error."

# --- CALLBACK FUNCTION FOR BUTTONS ---
def set_status(article_id, new_status):
    """Find the article in session_state and update its status."""
    for article in st.session_state.search_results:
        if article['id'] == article_id:
            article['status'] = new_status
            break

# --- USER INTERFACE ---
st.title("🔬 Manara AI - Screening Workspace")

# --- SIDEBAR FOR SEARCH AND SCOREBOARD ---
with st.sidebar:
    st.header("Search Controls")
    p = st.text_input("**P** - Population/Problem")
    i = st.text_input("**I** - Intervention")
    c = st.text_input("**C** - Comparison (Optional)")
    o = st.text_input("**O** - Outcome (Optional)")
    saudi_filter = st.checkbox("Filter for Saudi Arabia only", value=True)

    if st.button("Search PubMed", type="primary"):
        query = construct_pico_query(p, i, c, o, saudi_filter=saudi_filter)
        if query:
            with st.spinner("Searching..."):
                search_pubmed(query)
        else:
            st.warning("Please enter at least one PICO term.")
    
    st.divider()
    st.header("Screening Progress")
    # Calculate counts for the scoreboard
    if st.session_state.search_results:
        included_count = sum(1 for a in st.session_state.search_results if a['status'] == 'included')
        excluded_count = sum(1 for a in st.session_state.search_results if a['status'] == 'excluded')
        pending_count = sum(1 for a in st.session_state.search_results if a['status'] == 'pending')
        
        st.metric("✅ Included", included_count)
        st.metric("❌ Excluded", excluded_count)
        st.metric("🤔 Pending", pending_count)

# --- MAIN SCREENING AREA ---
st.header("Search Results")

if st.session_state.error:
    st.error(st.session_state.error)
elif not st.session_state.search_results:
    st.info("Enter a query in the sidebar and click 'Search PubMed' to begin.")
else:
    # Display each article with Include/Exclude buttons
    for index, article in enumerate(st.session_state.search_results):
        # Use an expander to neatly show the abstract
        with st.expander(f"{article['status'].upper()}: {article['title']}"):
            # Fetch and display abstract
            abstract = get_abstract(article['id'])
            st.markdown(abstract)
            
            st.divider()
            
            # Use columns for button layout
            col1, col2 = st.columns(2)
            
            with col1:
                st.button("✅ Include", key=f"include_{article['id']}", on_click=set_status, args=(article['id'], 'included'), use_container_width=True)
            
            with col2:
                st.button("❌ Exclude", key=f"exclude_{article['id']}", on_click=set_status, args=(article['id'], 'excluded'), use_container_width=True)