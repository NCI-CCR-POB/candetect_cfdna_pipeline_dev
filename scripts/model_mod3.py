import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold
import timeit

RANDOM_SEED = 0
N_PARAMS = 20  # number of values for C to test
N_SPLITS = 5  # number of validation folds (k), original value is 5
N_REPEATS = 10  # number of k-fold validation iterations ##original value is 10

PARAM_GRID = {
    "logisticregression__C": [10**x for x in np.arange(1, 5 + ((5-1)/N_PARAMS), (5-1)/N_PARAMS)], # C = 1/lambda
    "logisticregression__penalty": ['l1'],  # alpha = 1 is lasso
    "logisticregression__class_weight": ['balanced'],
    "logisticregression__solver": ['saga'],  # faster on larger datasets
    "logisticregression__random_state": [RANDOM_SEED],
    "logisticregression__max_iter": [10000],
}

RSKF = RepeatedStratifiedKFold(
    n_splits=N_SPLITS, n_repeats=N_REPEATS, random_state=RANDOM_SEED)

def cross_validation(X, y, pipe, out_dir, name, param_grid=PARAM_GRID, splitter=RSKF, verbose=3, n_jobs=-1):
    print("\n\n---------------- Cross Validation ----------------\n")
    start_time = timeit.default_timer()

    clf = GridSearchCV(pipe, param_grid, cv=splitter,
                       verbose=verbose, n_jobs=n_jobs)
    clf.fit(X, y.values.ravel())

    best_clf = clf.best_estimator_
    print("\nBEST ESTIMATOR:")
    print(best_clf)

    elapsed = timeit.default_timer() - start_time
    print(f"Elapsed time: {elapsed}s", flush=True)
    return best_clf


def get_coeffs_and_probs(X, y, clf, splitter=RSKF,
                         model_name="logisticregression"):
    print("\n\n---------------- Data Collection ----------------\n")

    n_repeats = splitter.n_repeats
    n_splits = splitter.cvargs["n_splits"]
    best_coeffs_df = pd.DataFrame()
    prob_df = pd.DataFrame(index=y.index,
                           columns=[
                               f"repeat_{str(i)}" for i in range(n_repeats)]
                           )

    iter_idx = 0  # [0,N_SPLITS*N_REPEATS)
    repeat_idx = 0  # [0,N_REPEATS)
    split_idx = 0  # [0,N_SPLITS) - N_SPLITS per repeat
    for train_index, test_index in splitter.split(X, y.values.ravel()):
        iter_str = f"iter_{iter_idx}"
        print(iter_str)
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        clf.fit(X_train, y_train.values.ravel())

        # save coeffs - reorder so that zscores come first
        coeffs = clf[model_name].coef_[0]  # a numpy array
        # coeffs = np.concatenate([coeffs[-len(arms):], coeffs[:-len(arms)]])
        coeff_df = pd.DataFrame.from_dict({iter_str: coeffs})
        best_coeffs_df = pd.concat([best_coeffs_df, coeff_df], axis=1)

        # save preds
        y_probs = clf.predict_proba(X_test)
        test_libs = y_test.index.to_list()
        repeat_str = f"repeat_{repeat_idx}"
        for i in range(len(test_libs)):
            prob_df[repeat_str][test_libs[i]] = y_probs[i][1]

        iter_idx += 1
        split_idx += 1
        if split_idx == n_splits:
            repeat_idx += 1
            split_idx = 0

    return best_coeffs_df, prob_df




# # report coefficients for all cv iterations
# pcs = ["pc" + str(i) for i in range(best_coeffs_df.shape[0] - len(arms))]
# best_coeffs_df.set_index(pd.Index(arms+pcs), inplace=True)
# print(best_coeffs_df)
# best_coeffs_df.to_csv(f"{OUT_DIR}/{COMPARISON_NAME}.coeffs.csv")

# # report probs and delfi_score for each cv repeat
# prob_df["delfi_score"] = prob_df.mean(axis=1)
# prob_df = pd.concat([y, prob_df], axis=1)
# print(prob_df)
# prob_df.to_csv(f"{OUT_DIR}/{COMPARISON_NAME}.probs.csv")

# # report aucs from each cv repeat and delfi_score
# aucs = {}
# cols = prob_df.columns.to_list()[1:]
# y_true = prob_df["Diagnosis"].to_list()
# for col in cols:
#     y_score = prob_df[col].to_list()
#     aucs[col] = [roc_auc_score(y_true, y_score)]
# aucs_df = pd.DataFrame.from_dict(
#     aucs, columns=[COMPARISON_NAME], orient="index")
# print(aucs_df)
# aucs_df.to_csv(f"{OUT_DIR}/{COMPARISON_NAME}.aucs.csv")

# print("\n\n---------------- Done! ----------------\n")
