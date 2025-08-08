"""
Simple HTML Reporter for edna2obis
Uses only built-in Python libraries for cross-platform compatibility
"""
import os
import datetime
import webbrowser
from pathlib import Path


class HTMLReporter:
    def __init__(self, filename="edna2obis_report.html"):
        self.filename = filename
        self.sections = []
        self.status = "RUNNING"
        self.start_time = datetime.datetime.now()
        self.error_message = None
        self.warnings = [] # To track warning messages
        
    def _get_status_color(self):
        if self.status == "SUCCESS":
            return "#28a745"  # Green
        elif self.status == "FAILED":
            return "#dc3545"  # Red
        elif self.status == "WARNING":
            return "#ffc107"  # Yellow for warning
        else:
            return "#6c757d"  # Grey for running or other states
    
    def add_section(self, title, level=2):
        """Add a section header"""
        self.sections.append({
            'type': 'section',
            'content': f"<h{level}>{title}</h{level}>"
        })
    
    def add_text(self, text):
        """Add plain text"""
        self.sections.append({
            'type': 'text',
            'content': f"<p>{text}</p>"
        })
    
    def add_list(self, items, title=None):
        """Add a list of items"""
        content = ""
        if title:
            content += f"<h4>{title}</h4>"
        content += "<ul>"
        for item in items:
            content += f"<li>{item}</li>"
        content += "</ul>"
        
        self.sections.append({
            'type': 'list',
            'content': content
        })
    
    def add_dataframe(self, df, title=None, max_rows=10):
        """Add a pandas DataFrame as HTML table"""
        content = ""
        if title:
            content += f"<h4>{title}</h4>"
        
        # Show shape info
        content += f"<p><strong>Shape:</strong> {df.shape[0]:,} rows × {df.shape[1]} columns</p>"
        
        # Convert DataFrame to HTML, limiting rows
        df_display = df.head(max_rows) if len(df) > max_rows else df
        table_html = df_display.to_html(classes="table table-striped", escape=False)
        
        if len(df) > max_rows:
            content += f"<p><em>Showing first {max_rows} rows of {len(df):,} total rows</em></p>"
        
        # Wrap table in scrollable container with scroll hint
        scroll_hint = '<p><small><em>💡 Tip: Scroll horizontally if table extends beyond view</em></small></p>' if df.shape[1] > 6 else ''
        content += f'<div class="table-container">{table_html}</div>{scroll_hint}'
        
        self.sections.append({
            'type': 'dataframe',
            'content': content
        })
    
    def add_success(self, message):
        """Add a success message"""
        self.sections.append({
            'type': 'success',
            'content': f'<div class="alert alert-success">{message}</div>'
        })
    
    def add_warning(self, message):
        """Add a warning message and track it"""
        self.warnings.append(message)
        self.sections.append({
            'type': 'warning',
            'content': f'<div class="alert alert-warning"><strong>WARNING:</strong> {message}</div>'
        })
    
    def add_error(self, message):
        """Add an error message and automatically set report status to FAILED"""
        self.error_message = message
        self.status = "FAILED"  # Automatically mark the report as failed
        self.sections.append({
            'type': 'error',
            'content': f'<div class="alert alert-danger"><strong>ERROR:</strong> {message}</div>'
        })
    
    def add_text_with_submission_logos(self, text):
        """Add text (logos are now in header, so this just adds styled text)"""
        content = f"""
        <p style="font-size: 18px; color: #495057; margin-bottom: 25px; font-weight: 500;">{text}</p>
        """
        self.sections.append({
            'type': 'text_with_logos',
            'content': content
        })
    
    def update_dataframe_from_file(self, title_identifier, filepath):
        """
        Find a previously added dataframe by its title and update it by re-reading from a file.
        This is crucial for refreshing data in the report after it's been modified on disk.
        """
        import pandas as pd
        
        found_and_updated = False
        for section in self.sections:
            # Check if the section is a dataframe and if its title contains the identifier
            if section['type'] == 'dataframe' and title_identifier in section.get('content', ''):
                
                print(f"Found dataframe section '{title_identifier}' in report. Re-reading from {filepath}...")
                
                # Re-read the dataframe from the provided filepath
                try:
                    df = pd.read_csv(filepath)
                    
                    # Re-create the HTML content for this section, just like in add_dataframe
                    # This assumes the title is in an h4 tag from the original add_dataframe call
                    title_html = section['content'].split("</h4>")[0] + "</h4>"
                    
                    new_content = title_html
                    new_content += f"<p><strong>Shape:</strong> {df.shape[0]:,} rows × {df.shape[1]} columns</p>"
                    table_html = df.to_html(classes="table table-striped", escape=False)
                    scroll_hint = '<p><small><em>💡 Tip: Scroll horizontally if table extends beyond view</em></small></p>' if df.shape[1] > 6 else ''
                    new_content += f'<div class="table-container">{table_html}</div>{scroll_hint}'
                    
                    # Update the section's content in place
                    section['content'] = new_content
                    print(f"Successfully updated report section for '{title_identifier}'.")
                    found_and_updated = True
                    break # Stop after finding and updating the first match
                
                except Exception as e:
                    # If re-reading fails, add a warning in the report instead of crashing
                    error_html = f'<div class="alert alert-danger"><strong>ERROR:</strong> Could not re-read and update report for {title_identifier}. Error: {e}</div>'
                    section['content'] += error_html
                    print(f"Failed to update report section: {e}")

        if not found_and_updated:
            print(f"Warning: Did not find a report section to update for '{title_identifier}'.")

    def set_success(self):
        """Mark the report as successful (only if not already failed)"""
        if self.status != "FAILED":
            self.status = "SUCCESS"
            
    def set_warning(self):
        """Mark the report with a warning status (only if not already failed)"""
        if self.status != "FAILED":
            self.status = "WARNING"
    
    def set_failed(self, error_message=None):
        """Mark the report as failed"""
        self.status = "FAILED"
        if error_message:
            self.error_message = error_message
    
    def set_status(self, status, error_message=None):
        """Set the status (SUCCESS or FAILED) with optional error message"""
        # Don't allow overriding FAILED status with SUCCESS or WARNING
        if self.status == "FAILED" and status in ["SUCCESS", "WARNING"]:
            return # Keep the FAILED status
            
        self.status = status
        if error_message:
            self.error_message = error_message
            self.add_error(error_message)
    
    def open_in_browser(self):
        """Open the HTML report in browser (writes file first)"""
        self._write_html()
        self._open_in_browser()
    
    def write_html(self):
        """Public method to write HTML (calls private method)"""
        self._write_html()
    
    def save_and_open(self):
        """Save the HTML file and open it in browser"""
        self._write_html()
        self._open_in_browser()
    
    def save(self):
        """Just save the HTML file without opening"""
        self._write_html()
    
    def _write_html(self):
        """Write the complete HTML file"""
        end_time = datetime.datetime.now()
        duration = end_time - self.start_time
        
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>edna2obis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f8f9fa; }}
        .container {{ max-width: 1200px; margin: 0 auto; background-color: white; padding: 40px; border-radius: 12px; box-shadow: 0 8px 16px rgba(0,0,0,0.1); border: 1px solid #e9ecef; }}
        .status {{ padding: 20px; margin-bottom: 20px; border-radius: 5px; text-align: center; font-size: 24px; font-weight: bold; color: white; background-color: {self._get_status_color()}; }}
        .alert {{ padding: 15px; margin: 10px 0; border-radius: 4px; }}
        .alert-success {{ background-color: #d4edda; border-color: #c3e6cb; color: #155724; }}
        .alert-warning {{ background-color: #fff3cd; border-color: #ffeaa7; color: #856404; }}
        .alert-danger {{ background-color: #f8d7da; border-color: #f5c6cb; color: #721c24; }}
        .table-container {{ overflow-x: auto; margin: 10px 0; max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }}
        .table-container::-webkit-scrollbar {{ height: 8px; }}
        .table-container::-webkit-scrollbar-track {{ background: #f1f1f1; border-radius: 4px; }}
        .table-container::-webkit-scrollbar-thumb {{ background: #c1c1c1; border-radius: 4px; }}
        .table-container::-webkit-scrollbar-thumb:hover {{ background: #a8a8a8; }}
        .table {{ width: 100%; border-collapse: collapse; margin: 0; min-width: 600px; }}
        .table th, .table td {{ padding: 8px 12px; text-align: left; border: 1px solid #ddd; white-space: nowrap; }}
        .table-striped tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .table th {{ background-color: #e9ecef; font-weight: bold; }}
        .metadata {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        h1 {{ color: #212529; font-size: 2.5rem; font-weight: 700; text-align: center; margin: 30px 0 20px 0; text-shadow: 0 1px 2px rgba(0,0,0,0.1); }}
        h2, h3, h4 {{ color: #343a40; }}
        pre {{ background-color: #f8f9fa; padding: 10px; border-radius: 4px; overflow-x: auto; }}
        .header-section {{ text-align: center; margin: 30px 0 50px 0; }}
        .logo-bar {{ display: flex; justify-content: center; align-items: center; gap: 40px; margin: 25px 0; opacity: 0.8; }}
        .logo-bar img {{ height: 60px; object-fit: contain; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="status">{self.status}</div>
        
        <div class="metadata">
            <h3>Run Information</h3>
            <p><strong>Start Time:</strong> {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>End Time:</strong> {end_time.strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Duration:</strong> {str(duration).split('.')[0]}</p>
            <p><strong>Report File:</strong> {self.filename}</p>
        </div>
        
        <div class="header-section">
            <h1>edna2obis Processing Report</h1>
            <div class="logo-bar">
                <img src="../images/aoml_logo.png" alt="AOML">
                <img src="../images/noaa_omics_logo.png" alt="NOAA Omics">
                <img src="../images/obis_logo.png" alt="OBIS">
                <img src="../images/gbif_logo.png" alt="GBIF">
            </div>
        </div>
"""
        
        # Add all sections
        for section in self.sections:
            html_content += section['content'] + "\n"
        
        html_content += """
    </div>
</body>
</html>
"""
        
        with open(self.filename, 'w', encoding='utf-8') as f:
            f.write(html_content)
    
    def _open_in_browser(self):
        """Open the HTML file in the default browser"""
        try:
            # Convert to absolute path for cross-platform compatibility
            file_path = Path(self.filename).resolve()
            webbrowser.open(f"file://{file_path}")
        except Exception as e:
            print(f"Could not open browser automatically: {e}")
            print(f"Please open {self.filename} manually in your browser")


# Simple test function
def test_reporter():
    """Test the HTML reporter with sample content"""
    reporter = HTMLReporter("test_report.html")
    
    reporter.add_section("Testing HTML Reporter", level=1)
    reporter.add_text("This is a test of the HTML reporting system.")
    
    # Test with a sample DataFrame
    import pandas as pd
    test_df = pd.DataFrame({
        'column1': ['A', 'B', 'C', 'D', 'E'],
        'column2': [1, 2, 3, 4, 5],
        'column3': ['test', 'data', 'for', 'html', 'display']
    })
    
    reporter.add_dataframe(test_df, "Sample DataFrame")
    
    # Test different message types
    reporter.add_success("This is a success message")
    reporter.add_warning("This is a warning message")
    
    # Test list
    reporter.add_list(['Item 1', 'Item 2', 'Item 3'], "Sample List")
    
    reporter.set_success()
    reporter.save_and_open()
    print("Test report created: test_report.html")


if __name__ == "__main__":
    test_reporter() 