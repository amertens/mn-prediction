"""Screenshots of the live dashboard for the deck (docs/slides/img/dashboard_*.png).

    python scripts/concept_slides/dashboard_screenshots.py

Needs playwright (python -m pip install playwright; python -m playwright install chromium).
Open the app a minute before running: the free shinyapps tier stops idle instances.
"""
import os
import time

from playwright.sync_api import sync_playwright

URL = "https://amertens.shinyapps.io/micronutrient-burden/"
OUT = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "slides", "img")


def main():
    with sync_playwright() as p:
        b = p.chromium.launch()
        pg = b.new_page(viewport={"width": 1600, "height": 900}, device_scale_factor=2)
        pg.goto(URL, wait_until="load", timeout=120000)
        pg.wait_for_selector("text=Where is deficiency?", timeout=120000); time.sleep(6)
        pg.screenshot(path=os.path.join(OUT, "dashboard_start.png"))
        pg.click("text=Where is deficiency?"); time.sleep(1)
        pg.click("text=Map explorer"); time.sleep(12)   # leaflet tiles and the district layer
        pg.screenshot(path=os.path.join(OUT, "dashboard_map.png"))
        b.close()
    print("wrote dashboard_start.png and dashboard_map.png in", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
