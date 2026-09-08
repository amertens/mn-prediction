"""One source of truth for the survey year of each country (see R/survey_years.R).

    from survey_years import SURVEY_YEAR          # {"Gambia": 2018, "Ghana": 2017, "Malawi": 2016, "SierraLeone": 2013}
    from survey_years import survey_years_table   # every row of metadata/survey_years.csv

The table is built from the interview dates of every dated cluster; survey_year
is the calendar year of the respondent-weighted median interview date.
"""
import csv
import os

_HERE = os.path.dirname(os.path.abspath(__file__))


def _find_table():
    d = _HERE
    for _ in range(6):
        f = os.path.join(d, "metadata", "survey_years.csv")
        if os.path.exists(f):
            return f
        d = os.path.dirname(d)
    raise FileNotFoundError("metadata/survey_years.csv not found above " + _HERE)


def survey_years_table():
    with open(_find_table(), encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def survey_years(protocol_only=True, lower=False):
    out = {}
    for r in survey_years_table():
        if protocol_only and str(r.get("in_protocol", "TRUE")).upper() != "TRUE":
            continue
        k = r["country"].lower() if lower else r["country"]
        out[k] = int(r["survey_year"])
    return out


SURVEY_YEAR = survey_years()
