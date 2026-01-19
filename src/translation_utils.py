"""
Translation utility for loading localized strings in Python modules.

This module provides functionality to load translation data from JSON files
and make it available to Python code that needs localized text strings.

Author: Alexandros Notas
Institution: National Technical University of Athens
"""

import json
import os

# Global translation objects
translations = {}
current_language = 'en'
available_languages = ['en', 'el']

def load_translations(language='en'):
    """
    Load translation data from JSON file.
    
    Args:
        language (str): Language code (default: 'en')
        
    Returns:
        dict: Translation data or empty dict if loading fails
    """
    global translations, current_language
    
    if language not in available_languages:
        print(f"Warning: Language {language} not supported. Using English fallback.")
        language = 'en'
    
    current_language = language
    
    try:
        # Get the directory of the current script
        current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        translation_file = os.path.join(current_dir, 'static', 'translations', f'{language}.json')
        
        with open(translation_file, 'r', encoding='utf-8') as f:
            translations[language] = json.load(f)
            
        return translations[language]
    except (FileNotFoundError, json.JSONDecodeError, KeyError) as e:
        print(f"Warning: Could not load translations for {language}: {e}")
        return {}

def load_all_translations():
    """
    Load all available translation files.
    
    Returns:
        dict: Dictionary with all loaded translations
    """
    global translations
    
    for language in available_languages:
        load_translations(language)
    
    return translations

def set_language(language):
    """
    Set the current language for translations.
    
    Args:
        language (str): Language code to set as current
        
    Returns:
        bool: True if language was set successfully, False otherwise
    """
    global current_language
    
    if language not in available_languages:
        print(f"Warning: Language {language} not supported.")
        return False
    
    if language not in translations:
        load_translations(language)
    
    current_language = language
    return True

def get_translation(key_path, fallback='', language=None):
    """
    Get a translation string using dot notation for nested keys.
    
    Args:
        key_path (str): Dot-separated path to translation key (e.g., 'thresholds.blast.glassWindows')
        fallback (str): Fallback text if translation not found
        language (str): Override language for this specific translation
        
    Returns:
        str: Translated text or fallback
    """
    target_language = language or current_language
    
    # Load translations if not already loaded
    if target_language not in translations:
        load_translations(target_language)
    
    def _resolve_translation(data, keys):
        current = data
        idx = 0
        while idx < len(keys):
            if not isinstance(current, dict):
                raise KeyError
            key = keys[idx]
            if key in current:
                current = current[key]
                idx += 1
                continue
            remaining = '.'.join(keys[idx:]).strip()
            if remaining and remaining in current:
                return current[remaining]
            raise KeyError
        return current

    keys = key_path.split('.')

    try:
        return _resolve_translation(translations[target_language], keys)
    except (KeyError, TypeError):
        # Try English fallback if current language fails
        if target_language != 'en':
            try:
                if 'en' not in translations:
                    load_translations('en')
                return _resolve_translation(translations['en'], keys)
            except (KeyError, TypeError):
                pass
        
        return fallback

def get_available_languages():
    """
    Get list of available language codes.
    
    Returns:
        list: List of available language codes
    """
    return available_languages.copy()

def get_current_language():
    """
    Get the current language code.
    
    Returns:
        str: Current language code
    """
    return current_language

# Load English translations on module import for compatibility
load_translations('en')
