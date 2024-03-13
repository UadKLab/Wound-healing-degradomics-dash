import pandas as pd
import numpy as np
import re 
from ast import literal_eval 
from dash import Dash, html, dcc, Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go


### Load data
protein_lookup_df = pd.read_csv('data/csv/Normalized_impaired_proteins.csv', delimiter=';', usecols=['Accession', 'Description'])
#peptide_lookup_df = pd.read_csv('data/csv/Normalized_impaired_peptides.csv', delimeter=';', usecols=['Annotated Sequence', 'Master Protein Accessions'])
protein_data_df = pd.read_csv('data/csv/4lowID_removed_joined_proteins.csv', delimiter=';', decimal=',')
#peptide_data_df = pd.read_csv('data/csv/4lowID_removed_joined_peptides.csv', delimeter=';', decimal=',')
patients_txt = open('data/patients.txt', 'r').read()


### Set-up
dtu_logo = 'assets/dtu_logo.png'

plot_margin = dict(
    l=80,
    r=80,
    b=20,
    t=20,  # top margin
    pad=8  # padding
)

# Patients groups (needs to be same as in patients.txt)
control = 'control'
impaired = 'impaired'
improving = 'improving'
chronic = 'chronic'

# Checklist for selecting patiens group/s (labels and values)
patients_checklist = [
    {'label': 'Control', 'value': control},
    {'label': 'Impaired', 'value': impaired},
    {'label': 'Improving', 'value': improving},
    {'label': 'Chronic', 'value': chronic}
]

# Line colors for patients groups
patients_colors = {
    control: '#4daf4a',
    impaired: '#984ea3',
    improving: '#377eb8',
    chronic: '#e41a1c'
}

### Data processing and transformation

## Utils

def group_data_df_transform_to_plot_format(group_data_df, group_column_name):
    '''
    Transforms a dataframe (protein_data_df or peptide_data_df) to be on a nice format to slice and plot.
    NaN values are ignored in aggregation.
    Use last experiment (T8) for all patient groups but set as NaN if no measurements.
    Values are transformed to log scale.
    Arguments: 
        group_data_df: pandas dataframe with raw input data
        group_column_name: string to identify the column name for the repsective group (Accession/Annotated Sequence)
    Returns:
        group_data_df: pandas dataframe with transformed input data
    '''
    # Set group as index and transpose
    group_data_df = group_data_df.set_index(group_column_name).T

    # Split index (time-patient_id) into time/experiment and patient_id multi-index
    group_data_df.index = group_data_df.index.str.split('-', expand=True)
    group_data_df.index.names = ['experiment', 'patient_id']

    # Map patient_id to groups
    group_data_df['patient_group'] = group_data_df.index.get_level_values('patient_id').map(lambda x: next((k for k, v in patients_ids.items() if x in v), None))

    # Group by time and patient_group
    group_data_df = group_data_df.groupby(['experiment', 'patient_group']).mean()

    # Reshape the df (unpivot group column) to be on a nice format to slice and plot the df 
    group_data_df = group_data_df.T.stack(level=0, future_stack=True).reset_index(inplace=False)
    group_data_df.columns.name = None

    # New column for the whole impaired patient group
    group_data_df[impaired] = group_data_df[[chronic, improving]].mean(axis=1)

    # Log transform
    group_data_df[[control,impaired,improving,chronic]] = group_data_df[[control,impaired,improving,chronic]].map(lambda x: np.log(x) if not np.isnan(x) and x > 0 else np.nan)

    # Confidence interval
    

    # Round to 1 decimal
    group_data_df = group_data_df.round(1)

    return group_data_df
    


## Dictionary with the patients groups as key and list of ids of patients as value
patients_ids = {}
patients_id_pattern = re.compile(r"(\w+)\s*=\s*\[([^\]]+)\]")   # remove dash in id if there to match with excel data
matches = patients_id_pattern.findall(patients_txt)
for key, value_str in matches:
    formatted_value_str = f"[{value_str.replace('-', '')}]"
    try:
        # Safely evaluate the string to a Python list
        value = literal_eval(formatted_value_str)
        patients_ids[key] = value
    except ValueError:
        print(f"Warning: Failed to parse '{formatted_value_str}' as a list.")


## Format peptide_lookup_df
#peptide_lookup_df['Master Protein Accessions'] = peptide_lookup_df['Master Protein Accessions'].str.split('; ')


## Transform data in protein_data_df and peptite_data_df - format so easy to plot and aggregate by time and patient groups
protein_data_df = group_data_df_transform_to_plot_format(protein_data_df, 'Accession')
#peptide_data_df = group_data_df_transform_to_plot_format(peptide_data_df, 'Annotated Sequence')


## 



## Clean-up after data processing
del patients_txt


# x-axis options
x_axis_options = protein_data_df['experiment'].unique().tolist()



#### Boiler-plate code ####

# Initialize Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
#load_figure_template('BOOTSTRAP')







# Define app layout
app.layout = html.Div(
    style={'width': '100%', 'height': '100vh'},
    children=[

        # Top banner with header and logo
        html.Div(
            style={'display': 'flex', 'justify-content': 'space-between', 
                   'align-items': 'center', 'text-align': 'center', 
                   'background-color':'#1c4c74', 'padding':'10px'},
            children=[
                html.Img(src=dtu_logo, style={'height':'55px'}),
                html.H1('Kostas Dashboard', style={'color':'white'})                
            ]
        ),

        # Search bar for protein and peptide
        html.Div(
            style={'margin-top':'20px', 'margin-bottom':'20px', 'margin-left': '5px'},
            children=[
                html.H3("Protein search"),
                html.Div(
                    style={'margin-left': '20px', 'margin-bottom':'10px', 'display': 'flex', 'align-items': 'center'},
                    children=[
                        html.Label("Description:" , style={'width':'100px'}),
                        dcc.Input(id='protein-description-input', 
                                  type='search', 
                                  debounce=True, 
                                  className='input-style'),
                    ]
                ),
                html.Div(
                    style={'margin-left': '20px', 'margin-bottom':'20px', 'display': 'flex', 'align-items': 'center'},
                    children=[
                        html.Label("Accession:", style={'width':'100px'}),
                        dcc.Dropdown(
                            className='dropdown-style',
                            id='protein-accession-dropdown',
                            options = [{'label': html.Div([html.B(row['Accession']), ' - ', row['Description']]), 'value': row['Accession']} for _, row in protein_lookup_df.iterrows()],
                            value=None
                        )
                    ]
                ),
                html.H3('Peptide search'),
                html.Div(
                    style={'margin-left': '20px', 'display': 'flex', 'align-items': 'center'},
                    children=[
                        html.Label("Sequence:" , style={'width':'100px'}),
                        dcc.Dropdown(id='peptide-sequence-dropdown',  
                                     className='dropdown-style'
                        ),
                    ]
                ),
            ]
        ),

        # Plot and filter by patients group/s
        html.Div(
            style={'flex':'1', 'display': 'flex', 'flex-direction': 'row', 
                   'align-items': 'center', 'margin-left': '5px'},
            children=[
                html.Div(
                    style={},
                    children=[
                        html.H5('Patient groups', style={'margin-bottom':'0px'}),
                        html.Div(
                            style={},
                            children=[
                                dcc.Checklist(
                                    id='patients-checklist-1',
                                    options=patients_checklist[:2],
                                    labelStyle={'display': 'block'}
                                )
                            ]
                        ),
                        html.Div(
                            style={'margin-left':'20px'},
                            children=[
                                dcc.Checklist(
                                    id='patients-checklist-2',
                                    options=patients_checklist[2:],
                                    labelStyle={'display': 'block'}
                                )
                            ]
                        )
                    ]
                ),
                html.Div(
                    style={'flex':'1'},
                    children=[
                        dcc.Tabs(
                            id='plot-tabs', 
                            value='protein-plot-tab',
                            className='plot-tabs-container',
                            children=[
                                dcc.Tab(
                                    label='Protein',
                                    value='protein-plot-tab',
                                    className='plot-tab',
                                ),
                                dcc.Tab(
                                    label='Peptide',
                                    value='peptide-plot-tab',
                                    className='plot-tab',
                                )
                            ]
                        ),
                        dcc.Graph(id='plot', style={'height':'60vh'})
                    ]
                )

                
            ]
        )
        
    ]
)

# TODO: gera peptide virkni
# TODO: setja restina í css










# Callback to update protein dropdown options based on search input
@app.callback(
    Output('protein-accession-dropdown', 'options'),
    [Input('protein-description-input', 'value')]
)
def update_dropdown(search_value):
    if search_value:
        # Filter options based on search value
        filtered_proteins = protein_lookup_df[protein_lookup_df['Description'].str.contains(search_value, case=False)]
        # Highlighting the matched parts
        highlighted_options = []
        for _, row in filtered_proteins.iterrows():
            highlighted = [html.B(row['Accession']), ' - ']
            segments = re.split(f'({search_value})', row['Description'], flags=re.IGNORECASE)
            for segment in segments:
                if re.search(search_value, segment, re.IGNORECASE):
                    highlighted.append(html.Span(segment, style={'background-color': 'yellow'}))
                else:
                    highlighted.append(html.Span(segment))
            highlighted_options.append({'label': html.Div(highlighted), 'value': row['Accession']})
        return highlighted_options
    else:
        return [{'label': html.Div([html.B(row['Accession']), ' - ', row['Description']]), 'value': row['Accession']} for _, row in protein_lookup_df.iterrows()]


# Callback to update protein dropdown options based on search input

# Callback to update peptide dropdown options based on search input
# @app.callback(
#     Output('peptide-dropdown', 'children'),
#     [Input('protein-dropdown', 'value'),
#      Input('peptide-description-input', 'value')]
# )
# def update_peptide_dropdown(selected_protein, search_term):
#     if not selected_protein:
#         return []

#     if search_term:
#         filtered_peptides = peptide_lookup_df[(peptide_lookup_df['Master Protein Accessions'].str.contains(selected_protein)) & 
#                                         (peptide_lookup_df['Annotated Sequence'].str.contains(search_term))]
#     else:
#         filtered_peptides = peptide_lookup_df[peptide_lookup_df['Master Protein Accessions'].str.contains(selected_protein)]

#     options = [{'label': row['Sequence'], 'value': row['Sequence']} for _, row in filtered_peptides.iterrows()]
#     return dcc.Dropdown(options=options)


# Callback to update plot based on selected protein and patients group/s
@app.callback(
    Output('plot', 'figure'),
    [Input('protein-accession-dropdown', 'value'),
     Input('patients-checklist-1', 'value'),
     Input('patients-checklist-2', 'value')]
)
def update_protein_plot(selected_protein, selected_groups_1, selected_groups_2):
    # Create the layout of the graph
        # Create the layout of the graph
    layout = go.Layout(
        title=None,
        plot_bgcolor='white',
        xaxis=dict(
            title='Experiment',
            gridcolor='lightgrey',
            tickmode='array',
            tickvals=list(range(len(x_axis_options))),
            ticktext=x_axis_options,
            categoryorder='array',
            categoryarray=x_axis_options,
            range=[-0.1, len(x_axis_options) - 0.9]  # Adjust the range to center the categories
        ),
        yaxis=dict(title='Log Transformed Values', type='linear', gridcolor='lightgrey'),
        showlegend=True,
        margin=plot_margin
    )


    if not selected_protein or (not selected_groups_1 and not selected_groups_2):
        traces = []
    
    else:
        # Filter
        filtered_protein = protein_data_df[protein_data_df['Accession'] == selected_protein]
        selected_groups = (selected_groups_1 or []) + (selected_groups_2 or [])

        # Traces
        traces = []
        for group in selected_groups:
            trace = go.Scatter(x=filtered_protein['experiment'], 
                            y=filtered_protein[group], 
                            mode='lines', 
                            name=group,
                            marker=dict(color=patients_colors[group]),
                            hovertemplate='Experiment: %{x}<br>Value: %{y}')
            traces.append(trace)
   
    # Create the figure
    fig = go.Figure(data=traces, layout=layout)

    return fig
    

if __name__ == '__main__':
    app.run_server(debug=True)