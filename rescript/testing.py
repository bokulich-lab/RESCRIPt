import requests


response = requests.get("https://www.bv-brc.org/api/genome_sequence/?in(genome_id,(224308.43))")

# Raise an error if the request was not successful
response.raise_for_status()

# Load the response data as JSON
data = response.json()

# Count the number of entries in the JSON dictionary
num_entries = len(data)

print(num_entries)