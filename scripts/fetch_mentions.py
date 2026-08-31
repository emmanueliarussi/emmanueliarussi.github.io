#!/usr/bin/env python3
"""Fetch a Google Alerts RSS/Atom feed and write data/mentions.json.

The feed URL is read from the GOOGLE_ALERTS_RSS_URL environment variable
(set as a repository secret and injected by the GitHub Action). If it is
missing, the existing data file is left untouched so the site keeps working.
"""

import html
import json
import os
import re
import sys
from datetime import datetime, timezone
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

FEED_URL = os.environ.get("GOOGLE_ALERTS_RSS_URL", "").strip()
OUT_PATH = os.path.join("data", "mentions.json")
MAX_ITEMS = 20
ATOM = "{http://www.w3.org/2005/Atom}"


def strip_html(text):
    return html.unescape(re.sub(r"<[^>]+>", "", text or "")).strip()


def real_url(google_redirect):
    """Google Alerts wraps links in a google.com/url?...&url=REAL redirect."""
    try:
        params = parse_qs(urlparse(google_redirect).query)
        if params.get("url"):
            return params["url"][0]
    except Exception:
        pass
    return google_redirect


def source_name(url):
    host = urlparse(url).netloc.lower()
    return host[4:] if host.startswith("www.") else host


def main():
    if not FEED_URL:
        print("GOOGLE_ALERTS_RSS_URL not set; leaving existing data unchanged.")
        return 0

    request = Request(FEED_URL, headers={"User-Agent": "Mozilla/5.0 (mentions-bot)"})
    with urlopen(request, timeout=30) as response:
        raw = response.read()

    root = ET.fromstring(raw)
    items = []
    for entry in root.findall(f"{ATOM}entry"):
        title = strip_html(entry.findtext(f"{ATOM}title", default=""))
        link_el = entry.find(f"{ATOM}link")
        href = link_el.get("href") if link_el is not None else ""
        url = real_url(href)
        published = entry.findtext(f"{ATOM}published", default="") or entry.findtext(
            f"{ATOM}updated", default=""
        )
        if not title or not url.lower().startswith(("http://", "https://")):
            continue
        items.append(
            {
                "title": title,
                "url": url,
                "source": source_name(url),
                "published": published,
            }
        )

    items = items[:MAX_ITEMS]
    os.makedirs("data", exist_ok=True)
    payload = {"updated": datetime.now(timezone.utc).isoformat(), "items": items}
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print(f"Wrote {len(items)} mentions to {OUT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
