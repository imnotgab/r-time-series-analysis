# PM2.5 Air Quality: Time Series Analysis & Forecasting
This repository contains an advanced time series analysis and forecasting project focused on PM2.5 particulate matter concentrations.

# Project Overview
The analysis explores how PM2.5 levels fluctuate over time, identifying underlying trends and seasonal patterns. Because environmental time series data often exhibits non-constant variance and strong seasonality, standard linear approaches can be insufficient. This project utilizes Box-Cox transformations and hybrid modeling techniques (combining regression with ARIMA errors) to uncover the true underlying dynamics and generate robust 24-month forecasts.

# Key Analytical Steps
## Data Preprocessing & Aggregation:

1. Converted raw daily PM2.5 measurements into monthly averages using the zoo package.

2. Split the dataset into a training set (up to June 2023) and a testing set (from July 2023 onwards) to validate model performance.

3. Exploratory Data Analysis (EDA):

4. Visualized seasonal patterns and trends using season plots, month plots, and lag plots.

5. Conducted Autocorrelation (ACF) and Partial Autocorrelation (PACF) analyses to identify seasonal lags and underlying autoregressive/moving average signatures.

6. Decomposed the time series to isolate trend, seasonal, and random components.

## Data Transformation & Stationarity:

1. Applied the Box-Cox transformation (with an optimized lambda) to the training set to stabilize variance and cleaned the data of extreme outliers.

2. Performed seasonal differencing (lag = 12) to achieve stationarity before ARIMA modeling.

## TSLM & Residual Modeling:

1. Built a baseline Time Series Linear Model (TSLM) capturing deterministic trend and seasonality.

2. Extracted TSLM residuals and modeled the remaining autocorrelation using explicit AR(12) and MA(12) processes.

3. Assessed coefficient significance and remodeled using fixed, statistically significant parameters.

## Advanced Forecasting Models:

1. Utilized the auto.arima algorithm on both the TSLM residuals and the base Box-Cox transformed data to find the optimal model order.

2. Applied Holt-Winters Exponential Smoothing as a robust alternative baseline.

3. Validated model residuals utilizing the Ljung-Box test and Shapiro-Wilk test to ensure white noise properties.

## Model Evaluation & Visualization:

1. Inverted the Box-Cox transformations for all forecasts to return them to their original scale.

2. Plotted all competing forecasts (TSLM+MA, TSLM+AR, Auto-ARIMA, Holt-Winters) simultaneously against the actual test set values.

3. Evaluated final model precision using the accuracy() function against the hold-out test set.

# Technologies Used
- Language: R

- Libraries: forecast, zoo

Remember to nstall any missing packages listed in the "Libraries" section.
