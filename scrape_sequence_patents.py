import requests
from bs4 import BeautifulSoup

def scrape_patent_by_id(patent_id):
    url = f"https://patents.google.com/patent/{patent_id}"

    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise HTTPError for bad responses (4xx or 5xx)
        soup = BeautifulSoup(response.content, 'html.parser')
        return soup
    except requests.exceptions.RequestException as e:
        print(f"Error fetching patent {patent_id}: {e}")
        return None


def extract_patent_data(patent_id):
    soup = scrape_patent_by_id(patent_id)
    if soup is None:
        return None

    # Title
    title_tag = soup.find("meta", {"name": "DC.title"})
    title = title_tag["content"] if title_tag else "Title not found"

    # Abstract
    abstract_tag = soup.find("meta", {"name": "DC.description"})
    abstract = abstract_tag["content"] if abstract_tag else "Abstract not found"

    # Extract claims from proper section
    claims = []
    claims_section = soup.find("section", itemprop="claims")
    if claims_section:
        claim_texts = claims_section.find_all("div", class_="claim-text")
        for claim_div in claim_texts:
            claims.append(claim_div.get_text(strip=True))

    extracted_info = {
        "patent_id": patent_id,
        "title": title,
        "abstract": abstract.replace('\n','').lstrip(),
        "claims": '\n'.join(claims)
    }
    return extracted_info


# Example usage
if __name__ == "__main__":
    patent_id = "US5817490"  # Replace with your patent ID
    data = extract_patent_data(patent_id)
    if data:
        for k,v in data.items():
            print(k)
            print(v)
            print()