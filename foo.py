import requests

r = requests.get("https://www.nrk.no")
print(r.status_code)
print(r.ok)
