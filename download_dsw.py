from datetime import datetime

def searchDSWx(ps):
    """
    Search for DSWx data based on given parameters.
    """
    # Check if date_start/date_end exist, otherwise extract from start/end parameters
    if not hasattr(ps, 'date_start') and hasattr(ps, 'start'):
        # Extract date from ISO format (2020-01-01T00:00:00Z) to YYYYMMDD
        try:
            # Remove quotes if present
            start_str = ps.start.strip("'\"")
            # Extract the date part (before T)
            date_part = start_str.split('T')[0]
            # Convert from YYYY-MM-DD to YYYYMMDD
            ps.date_start = date_part.replace('-', '')
            print(f"Extracted date_start {ps.date_start} from start parameter {ps.start}")
        except Exception as e:
            print(f"Error extracting date_start from start parameter: {e}")
            ps.date_start = '20200101'  # Default if extraction fails
            print(f"Using default date_start: {ps.date_start}")
    
    if not hasattr(ps, 'date_end') and hasattr(ps, 'end'):
        # Extract date from ISO format (2020-03-01T23:59:00Z) to YYYYMMDD
        try:
            # Remove quotes if present
            end_str = ps.end.strip("'\"")
            # Extract the date part (before T)
            date_part = end_str.split('T')[0]
            # Convert from YYYY-MM-DD to YYYYMMDD
            ps.date_end = date_part.replace('-', '')
            print(f"Extracted date_end {ps.date_end} from end parameter {ps.end}")
        except Exception as e:
            print(f"Error extracting date_end from end parameter: {e}")
            ps.date_end = '20200301'  # Default if extraction fails
            print(f"Using default date_end: {ps.date_end}")
    
    # If still missing, use defaults
    if not hasattr(ps, 'date_start'):
        ps.date_start = '20200101'
        print(f"No start date found. Using default date_start: {ps.date_start}")
    
    if not hasattr(ps, 'date_end'):
        ps.date_end = '20200301'
        print(f"No end date found. Using default date_end: {ps.date_end}")
    
    # Continue with existing code
    start_date = datetime(int(ps.date_start[0:4]), int(ps.date_start[4:6]), int(ps.date_start[6:8]))
    end_date = datetime(int(ps.date_end[0:4]), int(ps.date_end[4:6]), int(ps.date_end[6:8]))
    
    # Rest of the function remains unchanged 