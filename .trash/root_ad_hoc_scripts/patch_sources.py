import re

with open("src/screening_sources.py", "r") as f:
    code = f.read()

# Replace process_work_unit
old_func = """def process_work_unit(
    unit: WorkUnit,
    smiles_colname: str,
    chem_id_colname: str,
    archive_smiles_colname: str,
    archive_chem_id_colname: str,
) -> list[pd.DataFrame]:"""

new_func = """def process_work_unit(
    unit: WorkUnit,
    smiles_colname: str,
    chem_id_colname: str,
    archive_smiles_colname: str,
    archive_chem_id_colname: str,
    chunk_size: int = DEFAULT_INPUT_CHUNK_SIZE,
) -> Iterator[pd.DataFrame]:"""

code = code.replace(old_func, new_func)

code = re.sub(r'"""Parse a single WorkUnit into a normalized DataFrame."""\s+frames = \[\]', '"""Parse a single WorkUnit into a normalized DataFrame iterator."""', code)

impl = """
    if unit.source_type == "archive":
        # extract just this member from tar
        with tarfile.open(unit.source_path, "r:*") as archive:
            member = archive.getmember(unit.source_member)
            extracted = archive.extractfile(member)
            if extracted:
                row_offset = unit.start_row or 0
                for frame in pd.read_csv(TextIOWrapper(extracted, encoding="utf-8"), chunksize=chunk_size):
                    smiles_col = archive_smiles_colname if archive_smiles_colname in frame.columns else smiles_colname
                    chem_id_col = archive_chem_id_colname if archive_chem_id_colname in frame.columns else chem_id_colname
                    
                    if smiles_col not in frame.columns:
                        raise ValueError(f"Missing SMILES column in archive member {unit.source_member}")
                        
                    if chem_id_col not in frame.columns:
                        member_stem = _basename_without_suffixes(Path(unit.source_member))
                        frame[chem_id_col] = [f"{member_stem}__{row_offset + i + 1}" for i in range(len(frame))]
                        
                    frame = frame.rename(columns={smiles_col: "smiles", chem_id_col: "chem_id"})
                    frame["source_file"] = unit.source_member
                    frame["source_group"] = unit.source_group or unit.group_id
                    frame["source_member"] = unit.source_member
                    
                    if unit.start_row is not None and unit.max_rows is not None:
                        # Since we do not use skiprows for tar, slice manually or handle 
                        # Assuming start_row/max_rows isn't used for archive, we just assign source_row
                        pass
                    
                    frame["source_row"] = range(row_offset + 1, row_offset + 1 + len(frame))
                    row_offset += len(frame)
                    
                    frame = frame.dropna(subset=["smiles"]).copy()
                    if not frame.empty:
                        frame["smiles"] = frame["smiles"].astype(str)
                        frame["chem_id"] = frame["chem_id"].astype(str)
                        yield frame

    elif unit.source_type == "sqlite":
        with sqlite3.connect(unit.source_path) as conn:
            query = f'SELECT * FROM "{unit.source_member}"'
            if unit.start_row is not None and unit.max_rows is not None:
                query += f" LIMIT {unit.max_rows} OFFSET {unit.start_row}"
            
            row_offset = unit.start_row or 0
            for frame in pd.read_sql_query(query, conn, chunksize=chunk_size):
                if smiles_colname not in frame.columns:
                    raise ValueError(f"Missing SMILES column '{smiles_colname}' in {unit.source_path}")

                if chem_id_colname not in frame.columns:
                    frame[chem_id_colname] = [f"{unit.source_member}__{row_offset + 1 + i}" for i in range(len(frame))]

                frame = frame.rename(columns={smiles_colname: "smiles", chem_id_colname: "chem_id"})
                frame["source_file"] = str(unit.source_path)
                frame["source_group"] = unit.source_group or unit.group_id
                frame["source_member"] = unit.source_member
                frame["source_row"] = range(row_offset + 1, row_offset + 1 + len(frame))
                row_offset += len(frame)
                
                frame = frame.dropna(subset=["smiles"]).copy()
                if not frame.empty:
                    frame["smiles"] = frame["smiles"].astype(str)
                    frame["chem_id"] = frame["chem_id"].astype(str)
                    yield frame

    elif unit.source_type == "tabular":
        kwargs = {"sep": _infer_tabular_separator(Path(unit.source_path)), "chunksize": chunk_size}
        
        if unit.start_row is not None and unit.max_rows is not None and unit.start_row > 0:
            header_frame = pd.read_csv(unit.source_path, nrows=0, sep=kwargs["sep"])
            kwargs["skiprows"] = unit.start_row + 1
            kwargs["nrows"] = unit.max_rows
            kwargs["names"] = header_frame.columns
            kwargs["header"] = None
        elif unit.max_rows is not None:
            kwargs["nrows"] = unit.max_rows

        row_offset = unit.start_row or 0
        reader = pd.read_csv(unit.source_path, **kwargs)
        try:
            for frame in reader:
                if smiles_colname not in frame.columns:
                    raise ValueError(f"Missing SMILES column '{smiles_colname}' in {unit.source_path}")

                if chem_id_colname not in frame.columns:
                    frame[chem_id_colname] = [f"mol{row_offset + 1 + i}" for i in range(len(frame))]

                frame = frame.rename(columns={smiles_colname: "smiles", chem_id_colname: "chem_id"})
                frame["source_file"] = str(unit.source_path)
                frame["source_group"] = unit.source_group or unit.group_id
                frame["source_row"] = range(row_offset + 1, row_offset + 1 + len(frame))
                row_offset += len(frame)
                
                frame = frame.dropna(subset=["smiles"]).copy()
                if not frame.empty:
                    frame["smiles"] = frame["smiles"].astype(str)
                    frame["chem_id"] = frame["chem_id"].astype(str)
                    yield frame
        finally:
            reader.close()
"""

# Replace the body of process_work_unit
code = code[:code.find('    if unit.source_type == "archive":')] + impl

with open("src/screening_sources.py", "w") as f:
    f.write(code)
