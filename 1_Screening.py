# 1_Screening.py

import streamlit as st
import datetime
from Bio import Entrez
from sentence_transformers import SentenceTransformer, util
import re

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="Manara Screening", 
    page_icon="🔬", 
    layout="wide"
)

# --- NCBI API Configuration ---
Entrez.email = "your.email@example.com"

# --- MODEL LOADING ---
@st.cache_resource
def load_relevance_model():
    """Load the Sentence Transformer model for calculating relevance."""
    return SentenceTransformer("all-MiniLM-L6-v2")

model = load_relevance_model()

# --- INITIALIZE SESSION STATE ---
if "search_results" not in st.session_state:
    st.session_state.search_results = []
if "error" not in st.session_state:
    st.session_state.error = None
if "pico_query_text" not in st.session_state:
    st.session_state.pico_query_text = ""

# --- FUNCTIONS ---
def construct_pico_query(p=None, i=None, c=None, o=None, saudi_filter=False, year_range=None, publication_type=None):
    """Construct a high-precision PubMed query from all components."""
    terms = []
    if p: terms.append(f'"{p}"[Title/Abstract]')
    if i: terms.append(f'"{i}"[Title/Abstract]')
    if c: terms.append(f'"{c}"[Title/Abstract]')
    if o: terms.append(f'"{o}"[Title/Abstract]')

    if saudi_filter:
        terms.append(f'("Saudi Arabia"[MeSH Terms] OR "Saudi Arabia"[Title/Abstract])')

    if publication_type:
        publication_type_filter = [f'("{t}"[Publication Type])' for t in publication_type]
        publication_query = " OR ".join(publication_type_filter)
        terms.append(f'({publication_query})')

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
        abstract = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]
        return str(abstract)
    except (IndexError, KeyError, Exception):
        return None

def search_pubmed(query, pico_query_text, max_results=20): # Increased results for better screening
    """Search PubMed, calculate relevance, and update the session state."""
    st.session_state.pico_query_text = pico_query_text
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
        for record in records["PubmedArticle"]:
            article_id = record["MedlineCitation"]["PMID"]
            title = record["MedlineCitation"]["Article"]["ArticleTitle"]
            abstract = get_abstract(article_id)

            if abstract:
                results.append({
                    "title": title,
                    "id": article_id,
                    "status": "pending",
                    "abstract": abstract,
                    "checked": False # For multi-select
                })
                abstracts_to_score.append(abstract)

        if abstracts_to_score:
            query_embedding = model.encode(pico_query_text, convert_to_tensor=True)
            abstract_embeddings = model.encode(abstracts_to_score, convert_to_tensor=True)
            cosine_scores = util.cos_sim(query_embedding, abstract_embeddings)

            for i, article in enumerate(results):
                article["relevance_score"] = cosine_scores[0][i].item()

            results.sort(key=lambda x: x.get("relevance_score", 0), reverse=True)

        st.session_state.search_results = results
        st.session_state.error = None

    except Exception as e:
        st.session_state.search_results = []
        st.session_state.error = f"An API error occurred: {e}"

def highlight_text(text, keywords):
    """Highlights keywords in a case-insensitive manner in the given text."""
    if not isinstance(text, str) or not text.strip() or not isinstance(keywords, str) or not keywords.strip():
        return text

    keyword_list = [kw.strip() for kw in keywords.split() if kw.strip()]
    if not keyword_list:
        return text
        
    try:
        pattern = re.compile(r'\b(' + '|'.join(re.escape(kw) for kw in keyword_list) + r')\b', re.IGNORECASE)
        highlighted_text = pattern.sub(r'<mark>\1</mark>', text)
        return highlighted_text
    except re.error:
        return text

# --- CALLBACKS for multi-select ---
def update_checked_status():
    """Sync the 'checked' status in session_state from the checkbox widgets."""
    for article in st.session_state.search_results:
        widget_key = f"check_{article['id']}"
        if widget_key in st.session_state:
            article['checked'] = st.session_state[widget_key]

def bulk_set_status(new_status):
    """Set the status for all checked articles."""
    for article in st.session_state.search_results:
        if article.get('checked', False):
            article['status'] = new_status
            article['checked'] = False # Uncheck after action
            st.session_state[f"check_{article['id']}"] = False

def toggle_all_checkboxes():
    """Check or uncheck all article checkboxes based on the 'select all' widget."""
    is_checked = st.session_state.get("select_all_checkbox", False)
    for article in st.session_state.search_results:
        article['checked'] = is_checked
        st.session_state[f"check_{article['id']}"] = is_checked

def display_article(article):
    """Renders a single article card with a checkbox and expander."""
    score_percentage = article.get("relevance_score", 0) * 100
    
    col1, col2 = st.columns([0.05, 0.95])
    
    with col1:
        st.checkbox("select", key=f"check_{article['id']}", label_visibility="hidden")
    
    with col2:
        display_title = f"**{score_percentage:.1f}% Match** | {article['title']}"
        with st.expander(display_title):
            abstract = article.get("abstract", "Abstract not available.")
            highlighted_abstract = highlight_text(abstract, st.session_state.pico_query_text)
            st.markdown(highlighted_abstract, unsafe_allow_html=True)

# --- USER INTERFACE ---
st.title("🔬 Manara AI - Screening Workspace")

with st.sidebar:
    st.header("Search Controls")
    p = st.text_input("**P** - Population/Problem")
    i = st.text_input("**I** - Intervention")
    c = st.text_input("**C** - Comparison (Optional)")
    o = st.text_input("**O** - Outcome (Optional)")

    publication_type = st.multiselect(
        "Select the Publication Type",
        [
            "Randomized Controlled Trial",
            "Case Reports",
            "Case-control studies",
            "Cross-sectional studies",
            "Qualitative studies",
        ],
    )

    current_year = datetime.date.today().year
    year_range = st.slider(
        "Publication Year Range",
        min_value=1950,
        max_value=current_year,
        value=(current_year - 10, current_year),
    )

    saudi_filter = st.checkbox("Filter for Saudi Arabia only", value=True)

    if st.button("Search PubMed", type="primary"):
        pico_components = [p, i, c, o]
        pico_query_text = " ".join(filter(None, pico_components))

        query = construct_pico_query(
            p,
            i,
            c,
            o,
            saudi_filter=saudi_filter,
            year_range=year_range,
            publication_type=publication_type,
        )
        if query:
            with st.spinner("Searching PubMed & calculating relevance scores..."):
                # print(query)
                search_pubmed(query, pico_query_text)
        else:
            st.warning("Please enter at least one PICO term.")

    st.divider()
    st.header("Screening Progress")
    if st.session_state.search_results:
        included_count = sum(
            1 for a in st.session_state.search_results if a["status"] == "included"
        )
        excluded_count = sum(
            1 for a in st.session_state.search_results if a["status"] == "excluded"
        )
        # Calculate the 'maybe' count ---
        maybe_count = sum(
            1 for a in st.session_state.search_results if a["status"] == "maybe"
        )
        pending_count = sum(
            1 for a in st.session_state.search_results if a["status"] == "pending"
        )

        st.metric("✅ Included", included_count)
        st.metric("❌ Excluded", excluded_count)
        # Display the 'maybe' metric ---
        st.metric("🤔 Maybe", maybe_count)
        st.metric("⏳ Pending", pending_count)  # Changed emoji for clarity
        included_ids = [
            a["id"]
            for a in st.session_state.search_results
            if a["status"] == "included"
        ]
        if included_ids:
            st.success(
                f"Ready to move {len(included_ids)} articles to Data Extraction."
            )
            # This link will appear once you have at least one included article
            st.page_link(
                "pages/2_Data_Extraction.py", label="Go to Data Extraction", icon="➡️"
            )


st.header("Search Results")

if st.session_state.error:
    st.error(st.session_state.error)
elif not st.session_state.search_results:
    st.info("Enter your search criteria in the sidebar and click 'Search PubMed' to begin.")
else:
    update_checked_status() # Sync checkbox state at the start of the draw
    
    # --- Bulk Action Controls ---
    st.markdown("---")
    col1, col2, col3, col4, col5 = st.columns([2.5, 2.5, 2.5, 2, 3])
    with col1:
        st.button("✅ Include Selected", on_click=bulk_set_status, args=('included',), use_container_width=True)
    with col2:
        st.button("🤔 Maybe Selected", on_click=bulk_set_status, args=('maybe',), use_container_width=True)
    with col3:
        st.button("❌ Exclude Selected", on_click=bulk_set_status, args=('excluded',), use_container_width=True)
    with col5:
        st.checkbox("Select All / Deselect All", key="select_all_checkbox", on_change=toggle_all_checkboxes)
    st.markdown("---")
    
    # Filter articles into lists
    pending = [a for a in st.session_state.search_results if a['status'] == 'pending']
    included = [a for a in st.session_state.search_results if a['status'] == 'included']
    excluded = [a for a in st.session_state.search_results if a['status'] == 'excluded']
    maybe = [a for a in st.session_state.search_results if a['status'] == 'maybe']

    # Create the tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        f"⏳ Pending ({len(pending)})", 
        f"✅ Included ({len(included)})",
        f"❌ Excluded ({len(excluded)})",
        f"🤔 Maybe ({len(maybe)})"
    ])

    with tab1:
        st.subheader("Articles to be Screened")
        for article in pending:
            display_article(article)
    
    with tab2:
        st.subheader("Included Articles")
        for article in included:
            display_article(article)
            
    with tab3:
        st.subheader("Excluded Articles")
        for article in excluded:
            display_article(article)

    with tab4:
        st.subheader("Articles Marked as 'Maybe'")
        for article in maybe:
            display_article(article)