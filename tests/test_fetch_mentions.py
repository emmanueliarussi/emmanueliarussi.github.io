import io
from contextlib import redirect_stdout
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts import fetch_mentions as mentions


ATOM = b'''<feed xmlns="http://www.w3.org/2005/Atom"><entry>
  <title>Emmanuel &lt;b&gt;Iarussi&lt;/b&gt;</title>
  <link rel="self" href="https://example.org/feed-entry"/>
  <link href="https://www.google.com/url?url=https%3A%2F%2Fexample.org%2Farticle"/>
  <updated>2026-09-17T12:00:00Z</updated>
</entry></feed>'''
EMPTY = b'<feed xmlns="http://www.w3.org/2005/Atom"/>'


class FeedTests(unittest.TestCase):
    def setUp(self):
        # Expected failure cases must not create error annotations in Actions.
        self.enterContext(redirect_stdout(io.StringIO()))

    def test_atom_redirect_alternate_link_and_title(self):
        item, = mentions.parse_feed(ATOM)
        self.assertEqual(item, {
            "title": "Emmanuel Iarussi", "url": "https://example.org/article",
            "source": "example.org", "published": "2026-09-17T12:00:00Z",
        })

    def test_rss_dates_are_compatible_with_frontend(self):
        item, = mentions.parse_feed(b'''<rss><channel><item>
          <title>Interview</title><link>https://example.org/interview</link>
          <pubDate>Thu, 17 Sep 2026 12:00:00 GMT</pubDate>
        </item></channel></rss>''')
        self.assertEqual(item["published"], "2026-09-17T12:00:00+00:00")

    def test_reject_unexpected_xml_and_unusable_entries(self):
        for raw in (b"<html><body>Sign in</body></html>",
                    b'<feed xmlns="http://www.w3.org/2005/Atom"><entry/></feed>'):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                mentions.parse_feed(raw)

    def test_merge_preserves_history_deduplicates_and_sorts_before_limit(self):
        old = {"url": "https://example.org/old", "published": "2025", "title": "Old"}
        incoming = mentions.parse_feed(ATOM)
        self.assertEqual(mentions.merge_items([old], []), [old])
        self.assertEqual(len(mentions.merge_items([old] + incoming, incoming)), 2)
        with patch.object(mentions, "MAX_ITEMS", 1):
            self.assertEqual(mentions.merge_items(incoming, [old]), incoming)

    def test_empty_unchanged_invalid_and_missing_feed_preserve_file(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "mentions.json"
            original = json.dumps({"updated": "2026-09-01", "items": mentions.parse_feed(ATOM)})
            path.write_text(original)
            for raw, status in ((EMPTY, 0), (ATOM, 0), (b"<html/>", 1), (b"broken", 1)):
                with self.subTest(raw=raw), patch.object(mentions, "OUT_PATH", str(path)), \
                     patch.dict(os.environ, {"GOOGLE_ALERTS_RSS_URL": "https://example.org/feed"}), \
                     patch.object(mentions, "urlopen", return_value=io.BytesIO(raw)):
                    self.assertEqual(mentions.main(), status)
                    self.assertEqual(path.read_text(), original)
            with patch.dict(os.environ, {"GOOGLE_ALERTS_RSS_URL": ""}):
                self.assertEqual(mentions.main(), 1)
            self.assertEqual(path.read_text(), original)

    def test_new_mentions_are_written(self):
        with tempfile.TemporaryDirectory() as folder, \
             patch.object(mentions, "OUT_PATH", str(Path(folder) / "data/mentions.json")), \
             patch.dict(os.environ, {"GOOGLE_ALERTS_RSS_URL": "https://example.org/feed"}), \
             patch.object(mentions, "urlopen", return_value=io.BytesIO(ATOM)):
            self.assertEqual(mentions.main(), 0)
            self.assertEqual(len(json.loads(Path(mentions.OUT_PATH).read_text())["items"]), 1)


if __name__ == "__main__":
    unittest.main()
