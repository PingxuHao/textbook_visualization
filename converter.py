import markdown

# Read the content of the Markdown file
with open('chap_pde_1st.md', 'r') as md_file:
    md_content = md_file.read()

# Convert Markdown to HTML
html_output = markdown.markdown(md_content)

# Save the HTML output to a new file
with open('output.html', 'w') as html_file:
    html_file.write(html_output)
