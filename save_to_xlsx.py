import pandas as pd


def save_to_xlsx(xlsx_path, stresses_on_plane):
    # print(stresses_on_plane)

    data = []
    for plane in stresses_on_plane:
        data.append([
            plane.plane.dir, 
            plane.plane.dip, 
            plane.s_nn,
            plane.tau_n,
            ])
    df = pd.DataFrame.from_records(data)
    df.columns = ['dir', 'dip', 'σ_nn', '𝜏']
    df.to_excel(xlsx_path)
    # print(df)
