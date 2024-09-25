import pandas as pd
import numpy as np
from scipy import stats
import re 
from ast import literal_eval 
from dash import Dash, html, dcc, Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import os
from dotenv import load_dotenv
import boto3


### Environment variables
load_dotenv()
aws_access_key_id = os.getenv('AWS_ACCESS_KEY_ID')
aws_secret_access_key = os.getenv('AWS_SECRET_ACCESS_KEY')


### Load data
s3 = boto3.client(
    's3',
    aws_access_key_id=aws_access_key_id,
    aws_secret_access_key=aws_secret_access_key
)
bucket_name = 'wound-healing-storage'

obj = s3.get_object(Bucket=bucket_name, Key='Normalized_impaired_proteins.csv')
protein_lookup_df = pd.read_csv(obj['Body'], delimiter=';', usecols=['Accession', 'Description'])
obj = s3.get_object(Bucket=bucket_name, Key='Normalized_impaired_peptides.csv')
peptide_lookup_df = pd.read_csv(obj['Body'], delimiter=';', usecols=['Annotated Sequence', 'Modifications', 'Master Protein Accessions'])
obj = s3.get_object(Bucket=bucket_name, Key='Imputed_4lowID_removed_joined_proteins.csv')
protein_data_df = pd.read_csv(obj['Body'], delimiter=';', decimal=',')
obj = s3.get_object(Bucket=bucket_name, Key='Imputed_4lowID_removed_joined_peptides.csv')
peptide_data_df = pd.read_csv(obj['Body'], delimiter=';', decimal=',')
obj = s3.get_object(Bucket=bucket_name, Key='patients.txt')
patients_txt = obj['Body'].read().decode('utf-8')


### Set-up
dtu_logo = 'assets/dtu_logo.png'

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

# Colors for patient groups
patients_colors = {
    control: 'rgb(0.00784313725490196, 0.6196078431372549, 0.45098039215686275)',
    impaired: 'rgb(0.8, 0.47058823529411764, 0.7372549019607844)',
    improving: 'rgb(0.00392156862745098, 0.45098039215686275, 0.6980392156862745)',
    chronic: 'rgb(0.8352941176470589, 0.3686274509803922, 0.0)'
}
patients_colors_area = {
    control: 'rgba(0.00784313725490196, 0.6196078431372549, 0.45098039215686275, 0.20)',
    impaired: 'rgba(0.8, 0.47058823529411764, 0.7372549019607844, 0.20)',
    improving: 'rgba(0.00392156862745098, 0.45098039215686275, 0.6980392156862745, 0.20)',
    chronic: 'rgba(0.8352941176470589, 0.3686274509803922, 0.0, 0.20)'
}

plot_margin = dict(
    l=80,
    r=80,
    b=20,
    t=20,  # top margin
    pad=8  # padding
)

confidence = 0.75
confidence_print = int(100*confidence)
x_axis_options = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8']

### Data processing and transformation

## Utils

def confidence_interval(data, confidence):
    '''
    Calculates upper and lower confidence interval of series of data.
    Arguments:
        data: series/array of numeric data
        confidence: size of confidence interval [0,1]
    Returns:
        Upper and lower CI of data
    '''
    n = np.sum(~np.isnan(data))
    mean = np.nanmean(data)
    std_err = stats.sem(data, nan_policy='omit')
    if n <= 30:
        margin_error = std_err * stats.t.ppf((1 + confidence) / 2, n - 1)
    else:
         margin_error = std_err * stats.norm.ppf((1 + confidence) / 2)
        
    return mean - margin_error, mean + margin_error


def group_data_df_transform_to_plot_format(data_df, group_column_name):
    '''
    Transforms a dataframe (protein_data_df or peptide_data_df) to be on a nice format to slice and plot.
    NaN values are ignored in aggregation.
    Impaired group is a weighted average of chronic and improving.
    Values are transformed to log scale.
    Arguments: 
        group_data_df: pandas dataframe with raw input data
        group_column_name: string to identify the column name for the repsective group (Accession/Annotated Sequence)
    Returns:
        group_data_df: pandas dataframe with transformed input data
    '''
    # Set group as index and transpose
    data_df = data_df.set_index(group_column_name).T

    # Split index (time-patient_id) into time/experiment and patient_id multi-index
    data_df.index = data_df.index.str.split('-', expand=True)
    data_df.index.names = ['experiment', 'patient_id']

    # Map patient_id to groups
    data_df['patient_group'] = data_df.index.get_level_values('patient_id').map(lambda x: next((k for k, v in patients_ids.items() if x in v), None))

    # df to long format
    data_df = pd.melt(data_df.reset_index(), id_vars=['patient_group', 'experiment', 'patient_id'], var_name=group_column_name, value_name='value')

    # Add column for impaired patient group (mean of chronic and improving)
    impaired_df = data_df[data_df['patient_group'].isin(['chronic', 'improving'])].groupby([group_column_name,'experiment', 'patient_id'])['value'].mean().reset_index()
    impaired_df['patient_group'] = 'impaired'
    data_df = pd.concat([data_df, impaired_df], ignore_index=True)   # append the impaired values back into the original dataframe
    
    # Log transform (to make data normal distributed)
    data_df[['value']] = data_df[['value']].map(lambda x: np.log(x) if not np.isnan(x) and x > 0 else np.nan)

    # Group data to calcualte mean and CI for each Accession, experiment and patient group
    grouped = data_df.groupby(['experiment', 'patient_group', group_column_name])['value']
    grouped_df = pd.DataFrame(index=grouped.mean().index)
    # Mean
    grouped_df['mean'] = grouped.mean()
    # Confidence interval
    conf_interval = grouped.apply(confidence_interval, confidence=confidence)
    grouped_df['lower_ci'] = conf_interval.apply(lambda x: x[0])
    grouped_df['upper_ci'] = conf_interval.apply(lambda x: x[1])

    # Round to 2 decimals
    grouped_df = grouped_df.round(2)

    grouped_df.reset_index(inplace=True)

    return grouped_df
    

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

## Format and filter (only N-Term modifications) peptide_lookup_df
peptide_lookup_df['Master Protein Accessions'] = peptide_lookup_df['Master Protein Accessions'].str.split('; ')
peptide_lookup_df['Modifications'] = peptide_lookup_df['Modifications'].fillna('')
peptide_lookup_df = peptide_lookup_df[peptide_lookup_df['Modifications'].str.contains('.*\[N-term\].*', case=False)]

## Clean-up after data processing
del patients_txt


## Dash app

# Initialize Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Define app layout
app.layout = html.Div(
    className='dashboard-container',
    children=[

        # Top banner with header and logo
        html.Div(
            className='top-banner-container',
            children=[
                html.Img(src=dtu_logo, className='logo'),
                html.H1('Wound Healing Degradomics ', className='header')
            ]
        ),   # end top-banner-container

        # Search bar for protein and peptide
        html.Div(
            className='search-bar-container',
            children=[
                html.H3("Protein search"),
                html.Div(
                    className='search-bar-row',
                    children=[
                        html.Label("Description:" , className='search-bar-label'),
                        dcc.Input(id='protein-description-input', 
                                  type='search', 
                                  debounce=True, 
                                  className='dcc-input-style'),
                    ]
                ),
                html.Div(
                    className='search-bar-row',
                    children=[
                        html.Label("Accession:", className='search-bar-label'),
                        dcc.Dropdown(
                            className='dcc-dropdown-style',
                            id='protein-accession-dropdown',
                            options = [{'label': html.Div([html.B(row['Accession']), ' - ', row['Description']]), 'value': row['Accession']} for _, row in protein_lookup_df.iterrows()],
                            value=None
                        )
                    ]
                ),
                html.H3('Peptide search'),
                html.Div(
                    className='search-bar-row',
                    children=[
                        html.Label("Sequence:" , className='search-bar-label'),
                        dcc.Dropdown(id='peptide-sequence-dropdown',  
                                     className='dcc-dropdown-style'
                        ),
                    ]
                ),
            ]
        ),   # end search-bar-container

        # Plot and filter by patients group/s
        html.Div(
            className='plot-container',
            children=[
                dcc.Tabs(
                    id='plot-tabs', 
                    value='protein-plot-tab',
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
                html.Div(
                    className='plot-and-checklist-container',
                    children=[
                        html.Div(
                            style={},
                            children=[
                                html.H5('Patient groups'),
                                html.Div(
                                    children=[
                                        dcc.Checklist(
                                            id='patients-checklist-1',
                                            options=patients_checklist[:2],
                                            labelClassName='checklist-label',
                                            inputClassName='checklist-input'
                                        )
                                    ]
                                ),
                                html.Div(
                                    className='checklist-impaired',
                                    children=[
                                        dcc.Checklist(
                                            id='patients-checklist-2',
                                            options=patients_checklist[2:],
                                            labelClassName='checklist-label',
                                            inputClassName='checklist-input'
                                        )
                                    ]
                                )
                            ]
                        ),

                        dcc.Graph(id='plot', className='plot-graph', config={'displayModeBar': False})

                    ]   
                )   # end plot & group_checklist container
        
            ]   
        )   # end whole plot-container
    
    ]   
)   # end dashboard-container



## Callbacks

# Callback to update protein dropdown options based on search input
@app.callback(
    Output('protein-accession-dropdown', 'options'),
    [Input('protein-description-input', 'value')]
)
def update_protein_dropdown(search_value):
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


# Callback to show peptides in dropdown that are part of the selected protein
@app.callback(
        Output('peptide-sequence-dropdown', 'options'),
        [Input('protein-accession-dropdown', 'value')]
)
def update_peptide_dropdown(selected_protein):
    if selected_protein:
        filtered_peptides = peptide_lookup_df[peptide_lookup_df['Master Protein Accessions'].apply(lambda x: selected_protein in x)]
        return [{'label': row['Annotated Sequence'], 'value': row['Annotated Sequence']} for _, row in filtered_peptides.iterrows()]
    else:
        return []


# Callback to update plot based on selected protein/peptide and patients group/s and protein/peptide tab
@app.callback(
    Output('plot', 'figure'),
    [Input('protein-accession-dropdown', 'value'),
     Input('peptide-sequence-dropdown', 'value'),
     Input('patients-checklist-1', 'value'),
     Input('patients-checklist-2', 'value'),
     Input('plot-tabs', 'value')]
)
def update_plot(selected_protein, selected_peptide, selected_groups_1, selected_groups_2, selected_tab):

    traces = []
    annotations = [
        dict(
            x=0.5,
            y=0.5,
            text="No Data",
            font=dict(size=48, color="black"),
            showarrow=False,
            xref="paper",
            yref="paper"
        )
    ]

    # Plot protein or peptide
    selected = None
    if (selected_protein and selected_tab=='protein-plot-tab'):
        selected = selected_protein
        # Transform data on the fly based on filter
        group_column_name = 'Accession'
        data_df = protein_data_df
    elif (selected_peptide and selected_tab=='peptide-plot-tab'):
        selected = selected_peptide
        group_column_name = 'Annotated Sequence'
        data_df = peptide_data_df

    # Which patient groups to plot
    if selected and (selected_groups_1 or selected_groups_2):         

        # Transform data on the fly based on filter
        filtered_data_df = group_data_df_transform_to_plot_format(data_df[data_df[group_column_name] == selected], group_column_name)

        selected_groups = (selected_groups_1 or []) + (selected_groups_2 or [])            
        
        if not filtered_data_df.empty:
            annotations = None

            # Traces
            traces = []
            for group in selected_groups:
                # CI traces
                lower_ci_trace = go.Scatter(
                    x=x_axis_options,
                    y=filtered_data_df[filtered_data_df['patient_group']==group]['lower_ci'].values,
                    mode='lines',
                    line=dict(color='rgba(0,0,0,0)'),  # Set line color to transparent
                    fillcolor=patients_colors_area[group],
                    hoverinfo='none',
                    showlegend=False
                )
                traces.append(lower_ci_trace)
            
                upper_ci_trace = go.Scatter(
                    x=x_axis_options,
                    y=filtered_data_df[filtered_data_df['patient_group']==group]['upper_ci'].values,
                    mode='lines',
                    line=dict(color='rgba(0,0,0,0)'),  # Set line color to transparent
                    fill='tonexty',
                    fillcolor=patients_colors_area[group],
                    hoverinfo='none',
                    showlegend=False
                )
                traces.append(upper_ci_trace)

                # Mean trace
                customdata_values = np.column_stack((filtered_data_df[filtered_data_df['patient_group']==group]['lower_ci'].values,
                                                     filtered_data_df[filtered_data_df['patient_group']==group]['upper_ci'].values,
                                                     [confidence_print]* (8 if group==control else 7) ))
                trace_mean = go.Scatter(
                    x=x_axis_options, 
                    y=filtered_data_df[filtered_data_df['patient_group']==group]['mean'].values, 
                    mode='lines', 
                    name=group,
                    marker=dict(color=patients_colors[group]),
                    hovertemplate='Experiment: %{x}<br>Mean Value: %{y:.2f}<br>%{customdata[2]}% CI: [%{customdata[0]:.2f}, %{customdata[1]:.2f}]',
                    customdata=np.where(np.isnan(customdata_values), 'NaN', customdata_values)
                )
                traces.append(trace_mean)


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
        yaxis=dict(title='Log Transformed Values', type='linear', gridcolor='lightgrey', rangemode='tozero'),
        showlegend=True,
        margin=plot_margin,
        annotations=annotations
    )
            
    # Create the figure
    fig = go.Figure(data=traces, layout=layout)

    return fig
    

## Run app
if __name__ == '__main__':
    #app.run_server(debug=True)
    #app.run_server(debug=False)
    app.run_server(host='0.0.0.0', port=8000)




# TODO: fyrir cloud deployment nota S3 bucket (breyta innlestri gagnanna, kannski hreinsa sma til þvi 
# vid þurfum ekki allar þessar csv skrar) og nota WSGI production server til ad keyra appid (samt paeling
# ad skoda adeins betur þvi oft i docker eda configuration þa eru "run" leidbeiningar og þar er notad 
# WSGI server svo kannski er ekki naudsynlegt hér?? Skoda betur!)