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

