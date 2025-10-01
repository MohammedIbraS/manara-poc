
import requests
from transformers import pipeline

# Step 1: Fetch an article from the PubMed API (using a sample ID)
# We'll search for articles related to "metformin for diabetes".
# Note: It's good practice to add your email for the API.
print("1. Fetching article from PubMed...")
pubmed_api_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=31155736&retmode=json"
headers = {'User-Agent': 'AIReviewPOC/1.0 (mo7eh7@gmail.com)'}
response = requests.get(pubmed_api_url, headers=headers)
data = response.json()

# Extract the title and abstract (abstract is often in the 'fulljournalname' for summaries)
# A more robust solution would use the efetch API for full abstracts.
uid = '31155736'
title = data['result'][uid]['title']
# This is a simplified way to get some text; a full abstract pull is more complex.
journal_name = data['result'][uid]['fulljournalname']
article_text = f"Title: {title}. Journal: {journal_name}."
print(f"   - Fetched: {title[:50]}...")


# Step 2: Load a pre-trained "zero-shot classification" model from Hugging Face
print("\n2. Loading AI model from Hugging Face...")
# This model is great because you don't need to train it for your specific labels.
classifier = pipeline("zero-shot-classification", model="facebook/bart-large-mnli")
print("   - Model loaded successfully.")


# Step 3: Define labels and classify the abstract
print("\n3. Classifying the article text...")
candidate_labels = ["relevant to diabetes treatment","endometrial cancer" ,"cardiology research", "neurology study", "irrelevant"]
result = classifier(article_text, candidate_labels)
print("   - Classification complete.")


# Step 4: Print the result
print("\n4. Displaying results:")
print(f"\nArticle Text: \"{article_text}\"")
print("\nAI Analysis:")
# The results are sorted from highest to lowest score by default.
for label, score in zip(result['labels'], result['scores']):
    print(f"   - Label: {label:<30} | Confidence Score: {score:.2%}")