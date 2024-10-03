#!/home/ehsan/anaconda3/bin/python3

import argparse
import sqlite3
import csv
import os
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Import CSV data into an SQLite database.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input CSV file.')
    parser.add_argument('-db', '--database', required=True,
                        help='SQLite database file.')
    parser.add_argument('-n', '--table-name', required=True,
                        help='Name of the table in the database.')
    parser.add_argument('--force', action='store_true',
                        help='Force operation by overwriting incompatible tables.')
    return parser.parse_args()

def infer_column_types(sample_row):
    column_types = {}
    for column, value in sample_row.items():
        if value.isdigit():
            column_types[column] = 'INTEGER'
        else:
            try:
                float(value)
                column_types[column] = 'REAL'
            except ValueError:
                column_types[column] = 'TEXT'
    return column_types

def main():
    args = parse_arguments()

    input_csv = args.input
    db_file = args.database
    table_name = args.table_name
    force = args.force

    if not os.path.isfile(input_csv):
        print(f"Error: Input file '{input_csv}' does not exist.")
        sys.exit(1)

    # Read CSV file and get headers
    with open(input_csv, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        fieldnames = reader.fieldnames
        if not fieldnames:
            print("Error: CSV file must have a header row.")
            sys.exit(1)
        sample_row = next(reader, None)
        if not sample_row:
            print("Error: CSV file is empty.")
            sys.exit(1)
        column_types = infer_column_types(sample_row)
        csvfile.seek(0)  # Reset reader to the beginning of the file

    # Connect to SQLite database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Check if table exists
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table_name,))
    table_exists = cursor.fetchone()

    # Include 'ID' in expected columns
    expected_columns = ['ID'] + fieldnames

    if table_exists:
        # Verify schema compatibility
        cursor.execute(f"PRAGMA table_info('{table_name}')")
        existing_columns = [info[1] for info in cursor.fetchall()]
        if existing_columns != expected_columns:
            if force:
                print(f"Warning: Table '{table_name}' schema is incompatible. Recreating table due to --force.")
                cursor.execute(f'DROP TABLE "{table_name}"')
                table_exists = False  # Since we've dropped it
            else:
                print(f"Error: Table '{table_name}' exists with a different schema. Use --force to overwrite.")
                conn.close()
                sys.exit(1)
        else:
            print(f"Table '{table_name}' exists and is compatible. Appending data.")
    else:
        # Create new table with 'ID' as primary key
        columns = ['"ID" INTEGER PRIMARY KEY AUTOINCREMENT']
        columns += [f'"{name}" {column_types[name]}' for name in fieldnames]
        columns_sql = ', '.join(columns)
        cursor.execute(f'CREATE TABLE "{table_name}" ({columns_sql})')
        print(f"Table '{table_name}' created with 'ID' as primary key.")

    # Add UNIQUE constraint to 'Email' column if table is newly created
    if not table_exists:
        if 'Email' in fieldnames:
            cursor.execute(f'CREATE UNIQUE INDEX "idx_{table_name}_Email" ON "{table_name}" ("Email")')
            print(f"Unique index on 'Email' column created.")
        else:
            print("Warning: 'Email' column not found. Duplicate checks will not be performed.")

    # Fetch existing 'Email' and 'Professor' values to check for duplicates
    existing_emails = set()
    existing_professors = set()
    if table_exists:
        if 'Email' in fieldnames:
            cursor.execute(f'SELECT "Email" FROM "{table_name}"')
            existing_emails = set(row[0] for row in cursor.fetchall())
        if 'Professor' in fieldnames:
            cursor.execute(f'SELECT "Professor" FROM "{table_name}"')
            existing_professors = set(row[0] for row in cursor.fetchall())

    # Prepare insert statement (excluding 'ID' since it auto-increments)
    placeholders = ', '.join(['?'] * len(fieldnames))
    quoted_fieldnames = ', '.join([f'"{field}"' for field in fieldnames])
    insert_sql = f'INSERT INTO "{table_name}" ({quoted_fieldnames}) VALUES ({placeholders})'

    # Insert data
    with open(input_csv, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        rows_to_insert = []
        duplicate_emails = 0
        repeated_professors = set()
        for row in reader:
            row_values = []
            email = row.get('Email', '').strip()
            professor = row.get('Professor', '').strip()

            # Check for duplicate 'Email'
            if email and email in existing_emails:
                duplicate_emails += 1
                continue  # Skip inserting this row
            else:
                existing_emails.add(email)

            # Check for repeated 'Professor'
            if professor and professor in existing_professors:
                repeated_professors.add(professor)
            else:
                existing_professors.add(professor)

            for field in fieldnames:
                value = row[field]
                if column_types[field] == 'INTEGER':
                    value = int(value) if value else None
                elif column_types[field] == 'REAL':
                    value = float(value) if value else None
                else:
                    value = value.strip() if value else None
                row_values.append(value)
            rows_to_insert.append(tuple(row_values))

        if rows_to_insert:
            cursor.executemany(insert_sql, rows_to_insert)
            conn.commit()
            print(f"Inserted {len(rows_to_insert)} new records into '{table_name}'.")
        else:
            print("No new records to insert.")

        if duplicate_emails > 0:
            print(f"Skipped {duplicate_emails} records due to duplicate 'Email' values.")

        if repeated_professors:
            print("Warning: The following 'Professor' values are repeated:")
            for prof in repeated_professors:
                print(f" - {prof}")

    conn.close()

if __name__ == '__main__':
    main()
