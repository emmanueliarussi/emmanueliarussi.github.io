# emmanueliarussi.github.io

Google Alerts mentions appear under **In the media → Articles & press**.
The **News** section is edited manually in `index.html`.

The `Update media mentions` workflow reads the repository Actions secret
`GOOGLE_ALERTS_RSS_URL` twice daily. Use the RSS delivery link from Google Alerts,
not the alert management page or an email link.

To diagnose missing updates, manually run the workflow and inspect **Fetch mentions
from Google Alerts**. It reports the feed format, entry count, and usable mentions.
An empty feed produces a warning; a missing secret, unsupported response, or
unusable feed entries fails the run without replacing existing data. If it reports
zero entries, open the RSS link directly and check the alert's query and delivery
settings in Google Alerts. Scheduling a fetch cannot create entries absent from
the feed.

The script merges mentions by article URL, retains the newest 20, and changes
`data/mentions.json` only when the collected items change. An empty feed preserves
previous mentions. Run the feed regression tests with:

```sh
python3 -m unittest discover -s tests -v
```
