#!/usr/bin/env python3
"""Fetch a Google Alerts RSS/Atom feed and write data/mentions.json.

The feed URL is read from the GOOGLE_ALERTS_RSS_URL environment variable
(set as a repository secret and injected by the GitHub Action). If it is
missing or invalid, the run fails and the existing data file stays untouched.
Previously collected mentions are retained when entries leave the feed.
"""

import html
import json
import os
import re
import sys
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

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


def parse_feed(raw):
    root = ET.fromstring(raw)
    if root.tag == f"{ATOM}feed":
        entries = root.findall(f"{ATOM}entry")
        feed_format = "Atom"
    elif root.tag == "rss" and root.find("channel") is not None:
        entries = root.findall("channel/item")
        feed_format = "RSS"
    else:
        raise ValueError("Response is not an Atom or RSS feed. Check the Alerts RSS URL.")

    items = []
    for entry in entries:
        if feed_format == "Atom":
            title = strip_html(entry.findtext(f"{ATOM}title", default=""))
            href = next((link.get("href", "") for link in entry.findall(f"{ATOM}link")
                         if link.get("rel", "alternate") == "alternate"), "")
            published = entry.findtext(f"{ATOM}published", default="") or entry.findtext(
                f"{ATOM}updated", default=""
            )
        else:
            title = strip_html(entry.findtext("title", default=""))
            href = entry.findtext("link", default="")
            published = entry.findtext("pubDate", default="")
            if published:
                try:
                    published = parsedate_to_datetime(published).isoformat()
                except (ValueError, TypeError, OverflowError):
                    published = ""
        url = real_url(href)
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

    print(f"{feed_format} feed: {len(entries)} entries, {len(items)} usable mentions.")
    if entries and not items:
        raise ValueError("Feed has entries but none have a usable title and HTTP(S) link.")
    if not entries:
        print("::warning::Feed contains no entries. Check the alert query and RSS delivery in Google Alerts.")
    return items


def merge_items(existing, incoming):
    by_url = {item["url"]: item for item in existing}
    by_url.update({item["url"]: item for item in incoming})
    return sorted(by_url.values(), key=lambda item: item.get("published", ""), reverse=True)[:MAX_ITEMS]


def main():
    feed_url = os.environ.get("GOOGLE_ALERTS_RSS_URL", "").strip()
    if not feed_url:
        print("::error::GOOGLE_ALERTS_RSS_URL not set; existing data unchanged.")
        return 1

    try:
        request = Request(feed_url, headers={"User-Agent": "Mozilla/5.0 (mentions-bot)"})
        with urlopen(request, timeout=30) as response:
            incoming = parse_feed(response.read())
        path = Path(OUT_PATH)
        existing = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"items": []}
        items = merge_items(existing["items"], incoming)
        if items == existing["items"]:
            print(f"No changes; retaining {len(items)} mentions.")
            return 0
        payload = {"updated": datetime.now(timezone.utc).isoformat(), "items": items}
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        temporary.replace(path)
    except HTTPError as exc:
        print(f"::error::Feed request failed (HTTP {exc.code}); existing data unchanged.")
        return 1
    except (URLError, TimeoutError, OSError, ValueError, ET.ParseError, KeyError, TypeError) as exc:
        # Do not print exceptions that could include the private feed URL or response.
        print(f"::error::Unable to update mentions ({type(exc).__name__}). Check the RSS URL and feed format; existing data unchanged.")
        return 1
    print(f"Wrote {len(items)} mentions to {OUT_PATH} ({len(incoming)} from this feed).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
