def extract_features_from_text(sentence, prompt_template_path, token, model_name="RedHatAI/Llama-3.1-8B-Instruct",temp_val = 0, top_p = 1.0):
    """
    Sends a sentence to the LLM and extracts the output dict.
    """
    prompt = create_prompt(sentence, prompt_template_path)

    chat_response = client.completions.create(
        model=model_name,
        prompt=prompt,
        max_tokens=token,
        temperature=temp_val,
        top_p=top_p,
        #stop=["\n\n", "}"]
    )

    output_text = chat_response.choices[0].text.strip()
    #print('-----Hide this print after use----')
    #print(output_text)
    llm_dict = extract_dict(output_text)  # Assuming extract_dict function parses the LLM output JSON into dict
    
    if llm_dict is None:
        llm_dict = {}
        #print('-----Remove the print command after use------')
        #print(sentence)

    return llm_dict