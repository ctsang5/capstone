# Data Dictionary: options_data.csv

## Overview

This dataset contains equity and index options data. After cleaning, the dataset has **11,038 observations** and **7 variables**.

**Source:** Bloomberg Terminal from the University of San Francisco

**Date range:** 2/20/2026 - 8/21/2026

## Variables

| Variable | Type | Description | Example |
|---|---|---|---|
| symbol | string | Option contract identifier including ticker, expiry, strike, and type | AAPL US 02/20/26 C100 Equity |
| strike | float | Strike price of the option contract in USD | 100.0 |
| iv | float | Implied volatility of the option, expressed as a percentage | 48.21 |
| last_price | float | Last traded price of the option in USD | 164.07 |
| expiry | string | Expiration date of the option contract (M/D/YYYY) | 2/20/2026 |
| type | string | Option type | Call |
| underlying_price | float | Current price of the underlying asset in USD | 263.878 |

## Cleaning Steps

1. Trimmed dataset to first 7 columns
2. Dropped rows with missing values (14,802 to 11,038 rows)
3. Renamed columns for clarity (e.g., `equity` to `symbol`, `volatil` to `iv`)

## Notes

- All observations in the current dataset are calls; no puts are included
- Implied volatility values for deep in-the-money options appear unusually high (e.g., 4821 for AAPL), which may warrant further investigation
- Expiry is stored as a string, not a datetime object
