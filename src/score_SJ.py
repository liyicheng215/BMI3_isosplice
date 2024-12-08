import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier  # Import RandomForestClassifier
from sklearn.metrics import classification_report
from imblearn.over_sampling import SMOTE  # For oversampling minority class
import numpy as np
import argparse
from sklearn.model_selection import cross_val_score


def train_and_score_sj(csv_cluster, csv_label, threshold=None):
    df = pd.read_csv(csv_cluster)
    label_df = pd.read_csv(csv_label)
    label = label_df.iloc[:, 2]
    df['label'] = label

    sj_diff = (df.iloc[:, 4] - df.iloc[:, 3]).abs()
    df.insert(5, 'SJ_diff', sj_diff)
    # Separate the features and labels, set missing labels to -1 and split the labeled data
    feature_columns = df.columns[5:-1]
    label_column = df.columns[-1]
    df[label_column] = df[label_column].fillna(-1)
    labeled_df = df[df[label_column] != -1]
    X_labeled = labeled_df[feature_columns]
    y_labeled = labeled_df[label_column]
    
    # Standardize features and handle class imbalance
    scaler = StandardScaler()
    X_labeled_scaled = scaler.fit_transform(X_labeled)
    smote = SMOTE(random_state=42)
    X_resampled, y_resampled = smote.fit_resample(X_labeled_scaled, y_labeled)
    
    # Train a Random Forest model
    model = RandomForestClassifier(
    n_estimators=300,
    max_depth=30,
    min_samples_leaf=2, 
    min_samples_split=3,
    max_features='sqrt',
    oob_score=True
    )
    model.fit(X_resampled, y_resampled)
    
    # Predict probabilities score for each SJ
    X_all = df[feature_columns]
    X_all_scaled = scaler.transform(X_all)
    y_pred_prob = model.predict_proba(X_all_scaled)[:, 1]
    df['SJ_score'] = y_pred_prob
    output_file = 'scored_' + csv_cluster
    
    # Apply threshold if specified
    if threshold is not None:
        print(f"Applying threshold: {threshold}")
        df['Score_label'] = df['SJ_score'].apply(lambda x: 1 if x >= threshold else 0)

    df.to_csv(output_file, index=False)

    # Cross-validation on the resampled data
    cv_scores = cross_val_score(model, X_resampled, y_resampled, cv=5, scoring='accuracy')
    print(f"Cross-validation scores: {[f'{score:.6f}' for score in cv_scores]}")
    print(f"Average CV score: {np.mean(cv_scores):.6f}")
    print(f"OOB Score: {model.oob_score_:.6f}")

    # Print classification report for labeled data (model evaluation)
    y_labeled_pred = model.predict(X_resampled)
    print(classification_report(y_resampled, y_labeled_pred, digits=6))
    print(f"Scoring complete. Scores saved to {output_file}")
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train a random forest model and score SJ data.')
    parser.add_argument('csv_file', type=str, help='Path to the CSV file containing SJ cluster.')
    parser.add_argument('csv_label', type=str, help='Path to the CSV file containing SJ label.')
    args = parser.parse_args()
    train_and_score_sj(args.csv_file, args.csv_label)
