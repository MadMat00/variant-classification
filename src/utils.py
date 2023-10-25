def check_integrity(dataframe,log):
    df = dataframe.copy()
    for col in df.columns:
        assert not df[col].isnull().all()
    log.write_log("Dataframe integrity check passed", level="SUCCESS")