# pages/2_Data_Extraction.py

import streamlit as st
from Bio import Entrez
import requests
import base64
import fitz  # PyMuPDF
import json

# --- PAGE CONFIG ---
st.set_page_config(page_title="Data Extraction", page_icon="📝", layout="wide")

# --- API & URL Configuration (from secrets) ---
st.sidebar.header("Configuration")
try:
    # Priority 1: Try to load from secrets (most secure)
    n8n_webhook_url = st.secrets["N8N_WEBHOOK_URL"]
    st.sidebar.success("✅ n8n webhook configured securely.")
except (FileNotFoundError, KeyError):
    # Priority 2: Fallback to a user input if secrets are not found
    st.sidebar.warning("n8n webhook not found in secrets. Please provide it below.")
    n8n_webhook_url = st.sidebar.text_input(
        "Enter your n8n Webhook URL:",
        type="password", # Obscures the URL for security
        help="Paste the Production URL from your n8n webhook node here."
    )

# --- NCBI API Configuration ---
Entrez.email = "mo7eh7@gmail.com" 

# --- INITIALIZE SESSION STATE ---
if 'extraction_form_fields' not in st.session_state:
    st.session_state.extraction_form_fields = []
if 'current_article_id' not in st.session_state:
    st.session_state.current_article_id = None
if 'pdf_info' not in st.session_state:
    st.session_state.pdf_info = {}

# --- FUNCTIONS ---
def extract_text_from_pdf(pdf_source, is_url=True):
    """Extracts text from a PDF, either from a URL or an uploaded file object."""
    #st.write(pdf_source)
    try:
        if is_url:
            response = requests.get(pdf_source, timeout=15)
            response.raise_for_status()
            pdf_bytes = response.content
        else:
            pdf_source.seek(0)
            pdf_bytes = pdf_source.read()
        with fitz.open(stream=pdf_bytes, filetype="pdf") as doc:
            full_text = "".join(page.get_text() for page in doc)
        return full_text
    except Exception as e:
        return f"Error reading PDF: {e}"

def run_n8n_extraction(webhook_url, text_content, form_fields):
    """Sends data to an n8n webhook and robustly parses the AI response."""
    if not webhook_url:
        st.error("Please enter the n8n Webhook URL in your secrets file.")
        return {}

    headers = {"Content-Type": "application/json"}
    payload = {
        "text_content": text_content[:30000],
        "form_fields": form_fields
    }

    try:
        response = requests.post(webhook_url, headers=headers, json=payload, timeout=180)
        response.raise_for_status()
        
        response_data = response.json()
        
        # Navigate the complex structure from the n8n Gemini node
        text_from_ai = response_data['content']['parts'][0]['text']
        
        if not text_from_ai or not isinstance(text_from_ai, str):
            raise ValueError("The response from n8n was not the expected text format.")

        json_response = text_from_ai.strip().replace("```json", "").replace("```", "").strip()
        extracted_data = json.loads(json_response)
        return extracted_data
        
    except (requests.exceptions.RequestException, json.JSONDecodeError, KeyError, IndexError, ValueError) as e:
        st.error(f"An error occurred processing the n8n response: {e}")
        st.info("Ensure the n8n 'Respond to Webhook' node is configured correctly.")
        return {}
        
@st.cache_data
def get_doi_and_pdf_link(article_id):
    """Fetch an article's DOI, find its PDF link, and check if it's a download-only file."""
    error_message = None
    try:
        handle = Entrez.efetch(db="pubmed", id=article_id, rettype="fasta", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        doi = None
        for item in records['PubmedArticle'][0]['PubmedData']['ArticleIdList']:
            if item.attributes['IdType'] == 'doi':
                doi = str(item)
                break
    except Exception:
        return "Error fetching from PubMed", None, False, None

    if not doi:
        return "DOI not found", None, False, None

    unpaywall_url = f"https://api.unpaywall.org/v2/{doi}?email=mo7eh7@gmail.com"
    try:
        response = requests.get(unpaywall_url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        pdf_url = None
        best_location = data.get('best_oa_location')
        if best_location:
            pdf_url = best_location.get('url_for_pdf')
        
        if not pdf_url:
            oa_locations = data.get('oa_locations', [])
            for location in oa_locations:
                if location and location.get('url_for_pdf'):
                    pdf_url = location.get('url_for_pdf')
                    break 
        
        if not pdf_url:
            return doi, "No open-access PDF found", False, None

        is_attachment = False
        try:
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            with requests.get(pdf_url, stream=True, timeout=15, allow_redirects=True, headers=headers) as r:
                r.raise_for_status()
                if 'attachment' in r.headers.get('Content-Disposition', '').lower():
                    is_attachment = True
        except requests.exceptions.RequestException as e:
            error_message = f"Could not verify PDF type: {e}"
            is_attachment = True 

        return doi, pdf_url, is_attachment, error_message
            
    except requests.exceptions.RequestException:
        return doi, "Could not contact Unpaywall", False, None


# --- TEMPLATES ---
TEMPLATES = {
    "RCT": ["Sample Size", "Patient Population", "Intervention Group", "Control Group", "Primary Outcome"],
    "Observational": ["Sample Size", "Patient Population", "Exposure", "Outcome"],
    "Custom": []
}

# --- USER INTERFACE ---
st.title("📝 Data Extraction Workspace")

included_articles = [ a for a in st.session_state.get('search_results', []) if a.get('status') == 'included' ]

if not included_articles:
    st.warning("No articles have been marked as 'Included' in the Screening phase yet.")
    st.page_link("1_Screening.py", label="Back to Screening", icon="⬅️")
else:
    st.sidebar.header("Data Extraction Setup")
    template_choice = st.sidebar.selectbox("Choose a template:", options=list(TEMPLATES.keys()))

    if template_choice == "Custom":
        custom_fields_str = st.sidebar.text_area("Enter field names (one per line):", height=150)
        st.session_state.extraction_form_fields = [f.strip() for f in custom_fields_str.split('\n') if f.strip()]
    else: 
        st.session_state.extraction_form_fields = TEMPLATES[template_choice]

    st.sidebar.info(f"**Current Form Fields:** {', '.join(st.session_state.extraction_form_fields)}")

    st.info(f"You have moved {len(included_articles)} articles for data extraction.")
    st.markdown("---")
    
    selected_title = st.selectbox("Select an article:", options=[a['title'] for a in included_articles])
    selected_article = next((a for a in included_articles if a['title'] == selected_title), None)

    if selected_article:
        article_id = selected_article['id']
        
        if st.session_state.current_article_id != article_id:
            st.session_state.current_article_id = article_id
            with st.spinner("Searching for article PDF..."):
                doi, pdf_link, is_attachment, error_message = get_doi_and_pdf_link(article_id)
            st.session_state.pdf_info = {"doi": doi, "pdf_link": pdf_link, "is_attachment": is_attachment}
            #if error_message: st.toast(error_message, icon="⚠️")
        
        doi = st.session_state.pdf_info.get("doi")
        pdf_link = st.session_state.pdf_info.get("pdf_link")
        is_attachment = st.session_state.pdf_info.get("is_attachment")

        st.subheader(selected_title)
        st.caption(f"PubMed ID: {article_id} | DOI: {doi}")
        
        col1, col2 = st.columns(2)

        with col1:
            st.header("Article Viewer")
            viewable_pdf_url = pdf_link and "http" in pdf_link and not is_attachment
            uploaded_file = None
            if viewable_pdf_url:
                st.toast(" Open-access PDF found!", icon="✅")
                st.components.v1.iframe(pdf_link, height=800)
            else:
                if is_attachment:
                    st.toast(" PDF is a download-only file. Please use the download button or uploader.", icon="📄")
                    st.link_button("📥 Download PDF", url=pdf_link)
                if pdf_link == "No open-access PDF found" or pdf_link == "Could not contact Unpaywall":
                    st.warning("No open-access PDF found for this article.")
                    
                uploaded_file = st.file_uploader("Upload PDF:", type="pdf", key=f"uploader_{article_id}")
            
                if uploaded_file:
                    base64_pdf = base64.b64encode(uploaded_file.read()).decode('utf-8')
                    pdf_display = f'<iframe src="data:application/pdf;base64,{base64_pdf}" width="100%" height="800"></iframe>'
                    st.components.v1.html(pdf_display, height=800)

        with col2:
            st.header("Data Extraction Form")
            
            run_ai_button = st.button("🤖 Run AI Extraction via n8n", use_container_width=True)
            
            if run_ai_button:
                text_source, is_url_source = None, False
                if uploaded_file:
                    text_source, is_url_source = uploaded_file, False          
                elif viewable_pdf_url:
                    text_source, is_url_source = pdf_link, True
                
                if text_source and st.session_state.extraction_form_fields:
                    with st.spinner("Reading PDF and sending to n8n..."):
                        full_text = extract_text_from_pdf(text_source, is_url=is_url_source)
                        if "Error reading PDF" not in full_text:
                            extracted_data = run_n8n_extraction(n8n_webhook_url, full_text, st.session_state.extraction_form_fields)
                            
                            # Directly update the state of each widget
                            if extracted_data:
                                for field, value in extracted_data.items():
                                    widget_key = f"field_{article_id}_{field}"
                                    st.session_state[widget_key] = value
                        else: 
                            st.error(full_text)
                elif not st.session_state.extraction_form_fields:
                    st.warning("Please define fields in the sidebar.")
                else:
                    st.warning("A PDF must be available to run AI extraction.")

            if not st.session_state.extraction_form_fields:
                st.warning("Select a template or create a custom form in the sidebar.")
            else:
                for field in st.session_state.extraction_form_fields:
                    st.text_area(
                        field, 
                        key=f"field_{article_id}_{field}"
                    )
            
            if st.session_state.extraction_form_fields: 
                st.button("Save Data", type="primary")